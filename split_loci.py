import os
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Split .phy by RAD locus from VCF")
parser.add_argument("--vcf",  required=True, help="input VCF file")
parser.add_argument("--phy",  required=True, help="input .phy file")
parser.add_argument("-m", "--min-snp", type=int, default=1, help="minimum number of SNPs per locus (default: 1)")
parser.add_argument("-l", "--min-len", type=int, default=1, help="minimum sequence length per locus in bp (default: 1)")
parser.add_argument("-o", "--outdir",  default="loci_out", help="output folder name (default: loci_out)")
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

print("Loci passing filters: %d" % len(kept))

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

print("Loaded %d taxa x %d sites" % (len(taxa), len(seqs[0])))

# 4. write one .phy per locus
print("Writing loci to '%s/'..." % args.outdir)
for i, (lid, cols) in enumerate(sorted(kept.items(), key=lambda x: int(x[0])), 1):
    outpath = os.path.join(args.outdir, "locus_%s.phy" % lid)
    with open(outpath, "w") as out:
        out.write("%d %d\n" % (len(taxa), len(cols)))
        for name, seq in zip(taxa, seqs):
            sub = "".join(seq[c] for c in cols)
            out.write("%s  %s\n" % (name, sub))
    if i % 1000 == 0:
        print("  %d/%d done" % (i, len(kept)))

print("Done! %d loci saved to '%s/'" % (len(kept), args.outdir))
