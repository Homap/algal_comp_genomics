import argparse
from collections import Counter
from Bio import SeqIO
import statistics
from pathlib import Path

def calculate_n50(fasta_file):
    # 1. Get lengths
    lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    if not lengths:
        return 0
    # 2. Sort descending
    lengths.sort(reverse=True)

    # 3. Calculate total length and 50% threshold
    total_len = sum(lengths)
    half_sum = total_len / 2

    # 4. Find N50
    cumulative_sum = 0
    for length in lengths:
        cumulative_sum += length
        if cumulative_sum >= half_sum:
            return length

def fasta_stats(fasta_file):
    sequences = [str(record.seq).upper() for record in SeqIO.parse(fasta_file, "fasta")]
    # number of contigs/scaffolds/chromosomes
    num_contigs = len(sequences)
    # total assembly length
    assembly_len = sum(len(seq) for seq in sequences)
    # GC content
    base_count_list = [Counter(seq) for seq in sequences]
    total_base = sum(base_count_list, Counter())

    atgc_total = total_base["A"] + total_base["T"] + total_base["G"] + total_base["C"]
    GC_content = (total_base["G"] + total_base["C"])/ atgc_total if atgc_total > 0 else 0

    # per contig GC mean
    gc_per_contig = []
    for c in base_count_list:
        atgc = c["A"] + c["T"] + c["C"] + c["G"]
        if atgc > 0:
            gc_per_contig.append((c["G"] + c["C"]) / atgc)
    gc_mean = sum(gc_per_contig) / len(gc_per_contig) if gc_per_contig else 0
    gc_sd = statistics.stdev(gc_per_contig) if len(gc_per_contig) > 1 else 0

    # N_content
    total_bases = sum(total_base.values())
    N_content = (total_base["N"] / total_bases) if total_bases > 0 else 0

    return {
        "num_contigs": num_contigs,
        "assembly_len": assembly_len,
        "GC_content": GC_content,
        "gc_mean_per_contig": gc_mean,
        "gc_sd_per_contig": gc_sd,
        "N_content": N_content
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fasta_name = Path(args.fasta).name
    n50 = calculate_n50(args.fasta)
    stats = fasta_stats(args.fasta)

    header = [
        "assembly",
        "n50",
        "GC_content",
        "num_contigs",
        "assembly_len",
        "gc_mean_per_contig",
        "gc_sd_per_contig",
        "N_content"
    ]

    row = [
        fasta_name,
        n50,
        stats["GC_content"],
        stats["num_contigs"],
        stats["assembly_len"],
        stats["gc_mean_per_contig"],
        stats["gc_sd_per_contig"],
        stats["N_content"]
    ]

    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        out.write("\t".join(str(x) for x in row) + "\n")

if __name__ == "__main__":
    main()