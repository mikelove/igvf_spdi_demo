suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  cat("Script requires two arguments!\n")
  cat("Input variant file and name of output file\n")
  quit("no")
}
infile <- args[1]
outfile <- args[2]

################################ User input ###################################

max_len = 250 # maximum distance to search genome in overprecision algorithm

###############################################################################

# Read in list of variants
var <- read.table(infile, stringsAsFactors=F, header = F)
colnames(var) <- c("genome","ncbi_chr","hg38_pos","ref","alt")
var$seqnames <- sub("^0","",substr(var$ncbi_chr, 8,9))
var$seqnames <- sub("23","X",var$seqnames)
var$chr <- paste0("chr",var$seqnames)
var$ref_len <- nchar(var$ref) # Add length of ref in bp
var$alt_len <- nchar(var$alt) # Add length of alt in bp
var$indel_len <- var$ref_len - var$alt_len

# Get reference sequences around variants
X <- var %>%
  mutate(start = hg38_pos - max_len, end = hg38_pos + max_len)
X.gr <- makeGRangesFromDataFrame(X,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end")
X.seq <- getSeq(Hsapiens, X.gr)


# Function to identify insertion sequence
# in atypical ref/alt vcf format cases
get_indel = function(a, b) {
  a_len = nchar(a)
  b_len = nchar(b)
  
  if (a_len > b_len) {
    long = a
    long_end = a_len
    short = b
    short_end = b_len
  } else {
    long = b
    long_end = b_len
    short = a
    short_end = a_len
  }
  l = abs(a_len - b_len)
  start = 1
  
  while (substr(long, start, start) == substr(short, start, start)) {
    if (start == short_end) {
      return(list(substr(long, start+1, long_end), start+1))
    }
    start = start + 1
  }
  while (substr(long, long_end, long_end) ==  substr(short, short_end, short_end)) {
    long_end = long_end - 1
    short_end = short_end - 1
    if (long_end - start < l) {
      return(list(substr(long, start, long_end), start))
    }
  }
  if (long_end - start != l+1) {
    return(list("Error_Allele_Format", "Error_Allele_Format"))
  } else {
    return(list(substr(long, start, long_end), start))
  }
}


# Version of NCBI Variant Overprecision Correction Algorithm
speedy_name <- function(dnaseq, var) {
  seq_str <- as.character(dnaseq)
  ref <- substr(seq_str, max_len + 1, max_len + 1 + var$ref_len - 1)
  
  if (ref != var$ref) { 
    return(rep("Error_Ref_Mismatch", 2))
  } 
  
  if (var$ref_len == 1 & var$alt_len ==  1) {
    speedy <- paste(var$chr, var$hg38_pos, var$ref, var$alt, sep = ":")
    spdi <- paste(RefSeq_map[RefSeq_map$chr == var$chr, 2], 
      var$hg38_pos - 1, var$ref, var$alt, sep = ":")
    return(c(speedy, spdi))
  }

  if (var$ref_len == var$alt_len & var$ref_len > 1) {
    end <- nchar(var$ref)
    if (substr(var$ref, 1, 1) == substr(var$alt, 1, 1) |
        substr(var$ref, end, end) == substr(var$alt, end, end)) {
      return(rep("Error_Allele_Format", 2))
    } else {
      speedy <- paste(var$chr, var$hg38_pos, var$ref, var$alt, sep = ":")
      spdi <- paste(RefSeq_map[RefSeq_map$chr == var$chr, 2],
        var$hg38_pos - 1, var$ref, var$alt, sep = ":")
      return(c(speedy, spdi))
    }
  }
  
  l <- abs(var$indel_len)
  #indel <- max(var$ref, var$alt)
  #indel <- substr(indel, 2, l+1)
  
  if (var$ref_len > var$alt_len) {
    insert <- get_indel(var$ref, var$alt)
    indel <- insert[[1]]
    start <- insert[[2]]
    #indel <- substr(var$ref, 2, var$ref_len)
    k <- 1
    d <- max_len + start + l
    while (substr(indel, k, k) == substr(dnaseq, d, d)) {
      d <- d + 1
      k <- k + 1
      if (k == l+1) { k <- 1 }
    }
    k <- l
    u <- max_len + start
    while (substr(indel, k, k) == substr(dnaseq, u-1, u-1)) {
      u <- u - 1
      k <- k - 1
      if (k == 0) { k <- l }
    }
    speedy <- paste(var$chr, var$hg38_pos+u-(max_len + 1), 
                substr(dnaseq, u, d-1), substr(dnaseq, u, d-1-l),
                sep = ":")
    spdi <- paste(RefSeq_map[RefSeq_map$chr == var$chr, 2], 
              var$hg38_pos+u-(max_len + 1) - 1,
              substr(dnaseq, u, d-1), substr(dnaseq, u, d-1-l),
              sep = ":") 
    return(c(speedy, spdi))
  }
  
  if (var$ref_len < var$alt_len) {
    insert <- get_indel(var$ref, var$alt)
    indel <- insert[[1]]
    start <- insert[[2]]
    #indel <- substr(var$alt, 2, var$alt_len)
    subseq(dnaseq, start = max_len + 1, width = var$ref_len) <- var$alt
    k <- 1
    d <- max_len + start + l
    while (substr(indel, k, k) == substr(dnaseq, d, d)) {
      d <- d + 1
      k <- k + 1
      if (k == l+1) { k <- 1 }
    }
    k <- l
    u <- max_len + start
    while (substr(indel, k, k) == substr(dnaseq, u-1, u-1)) {
      u <- u - 1
      k <- k - 1
      if (k == 0) { k <- l }
    }
    speedy <- paste(var$chr, var$hg38_pos+u-(max_len + 1), 
                substr(dnaseq, u, d-1-l), substr(dnaseq, u, d-1),
                sep = ":")
    spdi <- paste(RefSeq_map[RefSeq_map$chr == var$chr, 2],
              var$hg38_pos+u-(max_len + 1) - 1,
              substr(dnaseq, u, d-1-l), substr(dnaseq, u, d-1),
              sep = ":")
    return(c(speedy, spdi))
  }
}


# Map between GRCh38 molecule name and RefSeq sequence name

RefSeq_start <- rep("NC_0000", 23)
RefSeq_mid <- c(paste0("0", as.character(1:9)),
  as.character(10:23))
RefSeq_end <- as.character(c(11, 12, 12, 12, 10, 12, 14, 11, 12, 11,
  10, 12, 11, 9, 10, 10, 11, 10, 10, 11, 9, 11, 11))
RefSeq <- paste0(RefSeq_start, RefSeq_mid, ".", RefSeq_end) 
RefSeq_map <- data.frame(chr = paste0("chr", c(1:22, "X")), 
                NC = RefSeq)

# Loop through variants in list
# This should be improved to avoid looping
# Current rate is ~ 40k variants / min.
n <- length(X.seq)
res <- var[,c("ncbi_chr","chr","hg38_pos","ref","alt")]
colnames(res) <- c("chrRefSeqID","chr","position","ReferenceAllele","AlternativeAllele")
res$SPDI <- NA
for (i in 1:n) {
  #res[i, 1] <- c(var[i, "input_variant"])
  res[i, "SPDI"] <- speedy_name(X.seq[i], var[i, ])[2]
}

write.table(res, file = outfile, row.names = F, col.names = T, quote = F, sep = "\t")
