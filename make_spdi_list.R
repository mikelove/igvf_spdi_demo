# this script takes as input a variant list that looks like "1_25253604_hg38_G_A"
# easily modifiable for other specifications, adjust `split_char` and `info_idx`
# note that this script assumes input is 1-based position
args <- commandArgs(trailingOnly=TRUE)

# what character separates the fields
split_char <- "_"
# which fields contain the sequence, position, ref, and alt
info_idx <- c(1, 2, 4, 5)

v <- scan(args[1], what="char", quiet=TRUE)
s <- strsplit(v, split=split_char)
cm <- read.delim("GRCh38p14_chrom_map.tsv")
chrom <- cm[,4]
names(chrom) <- cm[,1]

# SPDI uses RefSeq names for chromosomes, and is 0-based for position.
# assuming that input is 0-based we need to subtract 1
out <- sapply(s, function(x) paste(chrom[x[info_idx[1]]], # sequence
                                   as.numeric(x[info_idx[2]]) - 1, # position (0-based)
                                   x[info_idx[3]], # deleted
                                   x[info_idx[4]], # inserted
                                   sep=":"))

# output is 0-based and ready for NCBI Variant Services API
write(out, file="spdi_for_batch_processing.txt", ncolumns=1)
