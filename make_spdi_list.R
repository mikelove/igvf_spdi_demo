# this script takes as input a variant list that looks like "1_25253604_hg38_G_A"
# easily modifiable for other specifications, adjust `split_char` and `info_idx`
args <- commandArgs(trailingOnly=TRUE)

split_char <- "_"
info_idx <- c(1, 2, 4, 5)

v <- scan(args[1], what="char", quiet=TRUE)
s <- strsplit(v, split=split_char)
cm <- read.delim("GRCh38p14_chrom_map.tsv")
chrom <- cm[,4]
names(chrom) <- cm[,1]
# SPDI uses RefSeq names for chromosomes, and is 0-based for position
out <- sapply(s, function(x) paste(chrom[x[info_idx[1]]], # sequence
                                   as.numeric(x[info_idx[2]]) - 1, # position
                                   x[info_idx[3]], # deleted
                                   x[info_idx[4]], # inserted
                                   sep=":"))
write(out, file="spdi_for_batch_processing.txt", ncolumns=1)
