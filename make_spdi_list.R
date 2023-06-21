v <- scan("variants.txt", what="char", quiet=TRUE)
s <- strsplit(v, "_")
cm <- read.delim("GRCh38p14_chrom_map.tsv")
chrom <- cm[,4]
names(chrom) <- cm[,1]
# SPDI uses RefSeq names for chromosomes, and is 0-based for position
out <- sapply(s, function(x) paste(chrom[x[1]], as.numeric(x[2]) - 1,
                                   x[4], x[5], sep=":"))
write(out, file="spdi_for_batch_processing.txt", ncolumns=1)
