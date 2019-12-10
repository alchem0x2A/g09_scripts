#! /usr/bin/env python3
from src.gaussian import G09Calculator
from ase.parallel import parprint, paropen, rank, size, world
from argparse import ArgumentParser
import shutil, os
from pathlib import Path


def _parse_kw(string):
    try:
        key_dict = {ll[0]: ll[1] for ll in \
                    [l.strip().split("=") \
                     for  l in string.strip().split(",")]}
        return key_dict
    except (KeyError, AttributeError, IndexError):
        return dict()

def _cleanall(basedir):
    parprint("Clean up {0}".format(basedir))
    if rank == 0:
        p = Path(basedir)
        plist = []
        for type in ["*.com", "*.log", "*.ase",
                     "*.xyz", "*.chk", "lsf*"]:
            plist = plist + list(p.glob(type))

        for f in plist:
            if f.name != "init.xyz":
                os.unlink(f)
    

    
def main():
    parser = ArgumentParser()
    parser.add_argument("--process", "-p",
                        type=str,
                        default=None,
                        help="Process to be executed for a G09 calculation.")
    parser.add_argument("basedir",
                        type=str,
                        default=".",
                        help="Base directory for calculation. Should contain at least a init.xyz file")
    parser.add_argument("--config", "-c",
                        type=str,
                        help="Configuration json file containing basis set and defaults")
    parser.add_argument("--extras", "-e",
                        default="",
                        type=str,
                        help="Extra parameters separated by KEY1=VAL1,KEY2=VAL2")
    parser.add_argument("--clean", "-C",
                        default=False,
                        action='store_true',
                        help=('Clean up specific process. '
                              'If no process is specified, clean up .chk and .log files'))
    args = parser.parse_args()

    if args.process is None:
        if args.clean:
            _cleanall(args.basedir)
            return
        return
    g = G09Calculator(base=args.basedir,
                      config_file=args.config)
    g.process(label=args.process,
              clean=args.clean,
              **_parse_kw(args.extras))

if __name__ == "__main__":
    main()
