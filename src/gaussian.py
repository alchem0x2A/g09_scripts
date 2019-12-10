# Scripts for G09 calculations
from ase.io import read, write
from ase.parallel import parprint, paropen, rank, size, world
from ase.calculators.gaussian import Gaussian
from ase.io.gaussian import read_gaussian_out
from pathlib import Path
import json
from . import utils

proc_params = {"opt": dict(input="init.xyz",
                            output="opt.xyz",
                            kw=dict(opt="opt")),
                "gs": dict(input="opt.xyz",
                           inchk="opt",
                           output="gs.xyz",
                           kw=dict(population="reg",
                                   geom="allcheck",
                                   density="current",)),
                "td": dict(input="opt.xyz",
                           inchk="opt",
                           output="td-{root}.xyz",
                           kw=dict(td="(root={root},nstates={nstates})",
                                   opt="opt",
                                   density="current",
                                   population="reg",
                                   geom="allcheck"))}
class G09Calculator(Gaussian):
    """Modified class to wrap the original Gaussian class
       Each single calculation is done by separate method
    """
    def __init__(self, base=".",
                 mol="init.xyz",
                 config_file=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.load_config(config_file)
        self.base = Path(base)  # For later use
        try:
            mol_file = self.base / mol
            atoms = read(mol_file)
            self.init_atoms = atoms
        except Exception:
            raise

    def load_config(self, config_file=None):
        """`config_file contains the general params like Methods, basis-set, etc
        """
        if config_file is None:
            return              # Do nothing
        else:
            with open(config_file, "r") as f:
                params = json.load(f)
            self.set(**params)

    def _set_label(self, label=""):
        """Set the label under the base directory
        """
        l = self.base / label
        super().set_label(l.resolve().as_posix())

    def process(self, label, root=1, nstates=5, **kwargs):
        """Generalized process for optimization, 
           ground state, frequency and TD
        """
        if label not in proc_params.keys():
            raise KeyError("Process parameter not correct!")

        self._set_label(label)  # Causing writing label.com file
        params = proc_params[label]
        # in_atoms used to generate input
        try:
            in_file = self.base / params["input"]
            in_atoms = read(in_file)
        except Exception:
            raise

        out_files = [self.base / "{0}.chk".format(label),
                     self.base / "{0}.log".format(label)]
        if "output" in params.keys():
            out_struct = self.base / params["output"].format(root=root,
                                                             nstates=nstates)
            out_files.append(out_struct)
        else:
            out_struct = None
        # Calculation finished
        if all([o_.exists() for o_ in out_files]):
            parprint("Calculation process {0} is finished".format(label))
            #TODO: add post processing functions
            return True

        # Copy chk file from previous state
        if "inchk" in params.keys():
            utils.copy_chk(self.base,
                           params["inchk"], label)

        # Update parameters
        kw_mod = {k: v.format(root=root,
                              nstates=nstates,
                              **kwargs) for k, v in params["kw"].items()}
        self.set(**kw_mod)
        # Get environment
        nproc = utils.get_nproc()
        if nproc is None:
            nproc = 1
        mem = max(1024 * nproc, 8192)
        self.set(mem="{0}MB".format(mem),
                 nprocshared=nproc,
                 chk="{0}.chk".format(self.label))
        # Execute job now
        self.write_input(in_atoms)
        ec = utils.run_job(self.label)
        if ec != 0:             # Some error
            parprint("Process {0} failed, please check".format(label))
        out_atoms = read_gaussian_out(self.base / "{0}.log".format(label),
                                      quantity="structures")[-1]
        if out_struct is not None:
            write(out_struct, out_atoms)

        # return ec
