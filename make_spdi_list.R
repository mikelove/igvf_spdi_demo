v <- scan("variants.txt", what="char")
s <- strsplit(v, "_")
cm <- read.delim("GRCh38p14_chrom_map.tsv")
chrom <- cm[,4]
names(chrom) <- cm[,1]
out <- sapply(s, function(x) paste0(chrom[x[1]],":g.",x[2],x[4],">",x[5]))
write(out, file="spdi_for_batch_processing.txt", ncolumns=1)
