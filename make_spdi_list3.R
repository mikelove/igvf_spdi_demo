# this script takes as input a variant table with info in fields 2-5
# note that this script assumes input is 1-based position
args <- commandArgs(trailingOnly=TRUE)

# no scientific notation
options(scipen=999)

# which fields contain the sequence, position, ref, and alt
info_idx <- c(2, 3, 4, 5)

v <- read.table(args[1])

# SPDI uses RefSeq names for chromosomes, and is 0-based for position.
# assuming that input is 1-based we need to subtract 1
out <- paste(v[,info_idx[1]], # sequence
             as.numeric(v[,info_idx[2]]) - 1, # position (0-based)
             v[,info_idx[3]], # deleted
             v[,info_idx[4]], # inserted
             sep=":")

# output is 0-based and ready for NCBI Variant Services API
write(out, file="spdi_for_batch_processing.txt", ncolumns=1)
