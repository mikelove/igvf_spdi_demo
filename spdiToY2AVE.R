## SPDI to Y2AVE format converter
## Usage: Rscript spdiToY2AVE.R spdiFile 
##  spdiFile should have only 1 column (canonical SPDI), no header
##  Prints output table containing Y2AVE columns:
##    chrRefSeqID chr position ReferenceAllele AlternativeAllele SPDI
## Jesse Engreitz - Oct 9, 2023

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Script requires one argument: Input SPDI file\n")
  quit("no")
}
infile <- args[1]

# no scientific notation
options(scipen=999)

# Map between GRCh38 molecule name and RefSeq sequence name
RefSeq_start <- rep("NC_0000", 23)
RefSeq_mid <- c(paste0("0", as.character(1:9)),
  as.character(10:24))
RefSeq_end <- as.character(c(11, 12, 12, 12, 10, 12, 14, 11, 12, 11,
  10, 12, 11, 9, 10, 10, 11, 10, 10, 11, 9, 11, 11, 10))
RefSeq <- paste0(RefSeq_start, RefSeq_mid, ".", RefSeq_end) 
RefSeq_map <- data.frame(chr = paste0("chr", c(1:22, "X", "Y")), 
                NC = RefSeq)

if (FALSE) {

  ## Read input file
  spdi <- read.table(infile, stringsAsFactors=F)
  colnames(spdi) <- "SPDI"

  ## Process and format
  spdiSplit <- as.data.frame(do.call(rbind, strsplit(spdi$SPDI, ":")), stringsAsFactors=F)
  colnames(spdiSplit) <- c("chrRefSeqID","position","ReferenceAllele","AlternativeAllele")
  spdiSplit$SPDI <- spdi$SPDI
  spdiSplit <- merge(spdiSplit, RefSeq_map, by.x="chrRefSeqID", by.y="NC", all.x=TRUE)
  stopifnot(all(!is.na(spdiSplit$chr)))
  if (any(spdiSplit$ReferenceAllele == ""))
    spdiSplit$ReferenceAllele[spdiSplit$ReferenceAllele == ""] <- "-"
  if (any(spdiSplit$AlternativeAllele == ""))
    spdiSplit$AlternativeAllele[spdiSplit$AlternativeAllele == ""] <- "-"

  spdiSplit <- spdiSplit[,c("chrRefSeqID","chr","position","ReferenceAllele","AlternativeAllele","SPDI")]

  ## Write to stdout
  write.table(spdiSplit, sep='\t', quote=F, col.names=T, row.names=F)
}


# Open the input file
input_file <- file(infile, "r")

# Set the buffer size (N lines at a time)
buffer_size <- 100000  # Adjust as needed; buffered to avoid loading huge files into memory

# Read and process the file in chunks
while (TRUE) {
    lines_read <- readLines(input_file, n = buffer_size)
  
    # Check for end of file
    if (length(lines_read) == 0) {
      break
    }
  
    spdiSplit <- as.data.frame(do.call(rbind, strsplit(lines_read, ":")), stringsAsFactors=F)
    colnames(spdiSplit) <- c("chrRefSeqID","position","ReferenceAllele","AlternativeAllele")
    spdiSplit$SPDI <- lines_read
    spdiSplit <- merge(spdiSplit, RefSeq_map, by.x="chrRefSeqID", by.y="NC", all.x=TRUE)

    if (any(is.na(spdiSplit$chr))) {
      print("ERROR: Could not map chromosomes")
      print(subset(spdiSplit, is.na(chr)))
      quit("ERROR: Could not map chromosomes")
    }

    if (any(spdiSplit$ReferenceAllele == ""))
        spdiSplit$ReferenceAllele[spdiSplit$ReferenceAllele == ""] <- "-"
    if (any(spdiSplit$AlternativeAllele == ""))
        spdiSplit$AlternativeAllele[spdiSplit$AlternativeAllele == ""] <- "-"

    spdiSplit <- spdiSplit[,c("chrRefSeqID","chr","position","ReferenceAllele","AlternativeAllele","SPDI")]

    ## Write to stdout
    write.table(spdiSplit, sep='\t', quote=F, col.names=T, row.names=F)
}

# Close the input file
close(input_file)