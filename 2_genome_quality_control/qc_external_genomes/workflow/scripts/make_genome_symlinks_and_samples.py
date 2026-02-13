import os
from pathlib import Path
import argparse


# Author: Homa Papoli Yazdi
# This script produces symlink and sample table files for genomes downloaded from the NCBI.

def create_symlink(rootpath, linkpath, name_pattern, samples_tsv):
    genome_root = Path(rootpath)
    folders = [p for p in genome_root.iterdir() if p.is_dir()]

    link_dir = Path(linkpath)
    link_dir.mkdir(parents=True, exist_ok= True) # if parent directory do not exist, that is data/external_genomes, create them too.

    rows = []
    for f in folders:
        for i in f.glob(name_pattern):
            if i.name.startswith(f.name) and "cds_" not in i.name:
                dest = link_dir / f"{f.name}.fna"

                # remove existing link if present
                if dest.exists() or dest.is_symlink():
                    dest.unlink()

                # make the link target relative to the symlink directory
                rel_target = os.path.relpath(i, start=dest.parent) # os.path.relpath(i, start=dest.parent) computes the correct number of ../ hops from data/external_genomes/ to the real file
                dest.symlink_to(rel_target)

                # append to rows for the samples.tsv file
                rows.append((f.name, str(dest)))

    samples_tsv = Path(samples_tsv)
    samples_tsv.parent.mkdir(parents=True, exist_ok=True)

    with samples_tsv.open("w") as out:
        out.write("sample\tassembly\n")
        for sample, fasta in rows:
            out.write(f"{sample}\t{fasta}\n")
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