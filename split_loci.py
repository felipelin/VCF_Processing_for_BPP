import os
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Split .phy by RAD locus from VCF")
parser.add_argument("--vcf",     required=True, help="input VCF file")
parser.add_argument("--phy",     required=True, help="input .phy file")
parser.add_argument("-m", "--min-snp", type=int,   default=1,   help="minimum number of SNPs per locus (default: 1)")
parser.add_argument("-l", "--min-len", type=int,   default=1,   help="minimum locus length in bp (default: 1)")
parser.add_argument("-N", "--max-n",   type=float, default=0.0, help="maximum proportion of N per locus (0-1, default: 0 = no filter)")
parser.add_argument("-I", "--max-ind-n", type=float, default=0.0, help="maximum proportion of N per individual per locus (0-1, default: 0 = no filter)")
parser.add_argument("-o", "--outdir",  default="loci_out", help="output folder name (default: loci_out)")
parser.add_argument("--merge", action="store_true", help="also write all passing loci into one merged file")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# 1. read locus IDs and positions from VCF
print("Reading VCF locus IDs...")
locus_cols = defaultdict(list)
locus_pos  = defaultdict(list)
col = 0
with open(args.vcf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        locus_id = parts[2].split(":")[0]
        pos      = int(parts[1])
        locus_cols[locus_id].append(col)
        locus_pos[locus_id].append(pos)
        col += 1

# 2. filter by -m and -l
kept = {}
for lid in locus_cols:
    n_snp  = len(locus_cols[lid])
    length = max(locus_pos[lid]) - min(locus_pos[lid]) + 1
    if n_snp >= args.min_snp and length >= args.min_len:
        kept[lid] = locus_cols[lid]

print("  Total loci in VCF       : %d" % len(locus_cols))
print("  Loci passing -m and -l  : %d" % len(kept))

# 3. read .phy
print("Reading .phy file...")
taxa = []
seqs = []
with open(args.phy) as f:
    header = f.readline().split()
    n_taxa, n_sites = int(header[0]), int(header[1])
    for line in f:
        line = line.rstrip()
        if not line:
            continue
        parts = line.split(None, 1)
        taxa.append(parts[0])
        seqs.append(parts[1].replace(" ", ""))

print("  Loaded %d taxa x %d sites" % (len(taxa), len(seqs[0])))

# 4. apply -N filter and write loci
print("Writing loci to '%s/'..." % args.outdir)

n_filter = args.max_n
use_n_filter = (n_filter > 0.0)
n_removed_by_N = 0
i_filter = args.max_ind_n
use_i_filter = (i_filter > 0.0)
n_removed_by_I = 0
n_written = 0

merged_path = os.path.join(args.outdir, "merged.phy") if args.merge else None
merge_blocks = []  # list of (n_snp, list of seqs per taxon)

sorted_loci = sorted(kept.items(), key=lambda x: int(x[0]))

for i, (lid, cols) in enumerate(sorted_loci, 1):
    # extract sub-sequences for all taxa
    sub_seqs = ["".join(seq[c] for c in cols) for seq in seqs]

    # -N filter: check N proportion across all taxa x sites
    if use_n_filter:
        total_chars = sum(len(s) for s in sub_seqs)
        total_n     = sum(s.upper().count("N") for s in sub_seqs)
        n_prop = total_n / float(total_chars)
        if n_prop > n_filter:
            n_removed_by_N += 1
            continue

    # -I filter: check N proportion per individual
    if use_i_filter:
        failed = False
        for sub in sub_seqs:
            ind_n_prop = sub.upper().count("N") / float(len(sub))
            if ind_n_prop > i_filter:
                failed = True
                break
        if failed:
            n_removed_by_I += 1
            continue

    # write individual locus file
    outpath = os.path.join(args.outdir, "locus_%s.phy" % lid)
    with open(outpath, "w") as out:
        out.write("%d %d\n" % (len(taxa), len(cols)))
        for name, sub in zip(taxa, sub_seqs):
            out.write("^%s  %s\n" % (name, sub))

    if args.merge:
        merge_blocks.append((len(cols), sub_seqs))

    n_written += 1
    if i % 1000 == 0:
        print("  [%d/%d] processed..." % (i, len(kept)))

# report
print("")
print("=== Summary ===")
print("  Loci passing -m / -l    : %d" % len(kept))
if use_n_filter:
    print("  Removed by -N (>%.2f)   : %d" % (n_filter, n_removed_by_N))
if use_i_filter:
    print("  Removed by -I (>%.2f)   : %d" % (i_filter, n_removed_by_I))
print("  Loci written            : %d" % n_written)
print("  Output folder           : %s/" % args.outdir)

# 5. write merged file
if args.merge and merge_blocks:
    total_sites = sum(b[0] for b in merge_blocks)
    print("  Writing merged file     : %s" % merged_path)
    with open(merged_path, "w") as mf:
        for block_idx, (n_snp, block_seqs) in enumerate(merge_blocks):
            mf.write("%d %d\n" % (len(taxa), n_snp))
            for name, sub in zip(taxa, block_seqs):
                mf.write("^%-12s %s\n" % (name, sub))
            mf.write("\n")
    print("  Merged loci             : %d" % len(merge_blocks))
    print("  Merged total sites      : %d" % total_sites)

print("Done!")
