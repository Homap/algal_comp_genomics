import os
from pathlib import Path
import argparse

# Author: Homa Papoli Yazdi
# This script produces symlink and sample table files for genomes downloaded from the NCBI.

def create_symlink(rootpath, linkpath, name_pattern, samples_tsv):
    genome_root = Path(rootpath)

    # read sample names from samples_tsv (first column, skip header)
    with open(samples_tsv) as f:
        lines = f.readlines()[1:]
    samples = [line.strip().split()[0] for line in lines if line.strip()]

    link_dir = Path(linkpath)
    link_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for sample in samples:
        f = genome_root / sample   # folder corresponding to this sample

        for i in f.glob(name_pattern):
            if i.name.startswith(sample) and "cds_" not in i.name:
                dest = link_dir / f"{sample}.fna"

                if dest.exists() or dest.is_symlink():
                    dest.unlink()

                rel_target = os.path.relpath(i, start=dest.parent)
                dest.symlink_to(rel_target)

                rows.append((sample, str(dest)))

    return rows
"""
This pattern is about separating script mode from library mode.
If you run the file directly, it parses command-line args and runs. 
If Snakemake or another script imports your function, it does not run automatically.
"""

def main():
    ap = argparse.ArgumentParser() # creates a parser that understands command line arguments
    ap.add_argument("--rootpath", required=True)
    ap.add_argument("--linkpath", required=True)
    ap.add_argument("--pattern", required=True)
    ap.add_argument("--samples", required=True)
    args = ap.parse_args() # read the command line and store the values in args

    create_symlink(args.rootpath, args.linkpath, args.pattern, args.samples) # only run main() if this file is executed directly.

if __name__ == "__main__":
    main()