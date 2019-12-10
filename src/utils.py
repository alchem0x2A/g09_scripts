# Supporting scripts for job inputing stuffs
import os
import re
from ase.units import Hartree
import numpy as np
from pathlib import Path


def get_nproc():
    """Get max number of processes. If not in batch mode, return None
    """
    keyword = "LSB_MAX_NUM_PROCESSORS"
    try:
        n = os.environ[keyword]
        return int(n)
    except KeyError:
        return None


def confirm_scratch():
    """Confirm that `$TMPDIR is set correctly
    """
    if "TMPDIR" not in os.environ.keys():
        print("You are running Gaussian09 on Log-in node, this is not recommended.")
    else:
        tmpdir = os.environ["TMPDIR"]
        print("You are running Gaussian09 on Computing node. Will set the scratch.")
        os.environ["GAUSS_SCRDIR"] = tmpdir


def convert_cube(base, label, orbital="homo"):
    base = Path(base).resolve()
    # Convert chk file
    chkfile = base / "{0}.chk".format(label)
    fchkfile = base / "{0}.fchk".format(label)
    # cubfile, no syntax checking etc...
    cubfile = base / "{0}.cube".format(orbital)
    if not chkfile.exists():
        raise FileNotFoundError("No chk file, possibly not converged!")
    # Convert fchk
    if not fchkfile.exists():
        ec = os.system("formchk {0} {1}".format(chkfile, fchkfile))
        if ec != 0:
            raise ValueError("Error converting chk file to fchk!")
        ec1 = os.system("cubegen 1 MO=Homo {0} {1} -2".format(fchkfile, cubfile))  # fine grid
        if ec1 != 0:
            raise ValueError("Cubegen failed...")
    return True


def get_eigen(base, label):
    """Analysis eigenvalues from log file
    """
    base = Path(base)
    log_file = base / "{0}.log".format(label)
    if not log_file.exists():
        raise FileNotFoundError("No log file!")
    with open(log_file, "r") as f:
        text = f.read()
    # Population section
    pat_pop = r"\*{30,}[\s]+Population analysis.+\*{30,}(.+)Molecular Orbital Coeff"
    match = re.findall(pat_pop, text, re.DOTALL)  # match multiline
    if len(match) == 0:
        raise ValueError("Population calculation may be bad!")
    rest_text = match[0]
    pat_occ = r"Alpha\s+occ.\s+eigenvalues -- (.+)$"
    pat_unocc = r"Alpha\s+virt.\s+eigenvalues -- (.+)$"
    match_occ = re.findall(pat_occ, rest_text, re.MULTILINE)
    match_unocc = re.findall(pat_unocc, rest_text, re.MULTILINE)
    if (len(match_occ) == 0) or (len(match_unocc) == 0):
        raise ValueError("Matching eigenvalues failed!")
    eigen_occ = np.genfromtxt(["".join(match_occ)], delimiter=10)
    eigen_unocc = np.genfromtxt(["".join(match_unocc)], delimiter=10)
    homo = eigen_occ[-1] * Hartree
    lumo = eigen_unocc[0] * Hartree
    return homo, lumo

