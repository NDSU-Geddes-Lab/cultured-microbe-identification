#!/usr/bin/env Rscript

# 1. Create output table from input table with:
#      - original seq
#      - f_barcode, r_barcode, and well from matching seq against barcode table
#      - trimmed_seq
#      - one column for each plate+well combo, initialized with all zeros

# For each sequence in the seqtab:
#   1. Trim the primers
#   2. Identify the forward and reverse barcodes
#   3. Find the corresponding well number in the barcode table
#   4. Report the result in the column corresponding to the plate+well combo

seqtab <- read.csv("output/input_seqtab.csv")
output_exp <- read.csv("output/processed_data_output.csv")
barcodes <- read.csv("BC_to_well2.csv")

seq <- seqtab$X[1]

f_primer <- "GTGCCAGCMGCCGCGGTAA"
r_primer <- "GACTACHVGGGTATCTAATCC"

seq_len <- nchar(seq)

f_len <- nchar(f_primer)
r_len <- nchar(r_primer)

find_well <- function(seq) {
  i <- which(startsWith(seq, barcodes$F) & endsWith(seq, barcodes$R))
  
  if (length(i) != 1) {
    return(NA)
  } else {
    return(barcodes[i,3])
  }
}

wells <- unlist(lapply(seqtab$X, find_well)) # working

df <- data.frame(A=numeric(10), B=numeric(10), C=numeric(10))
df[,c("D","E","F")] <- 0

cols <- c("G","H","I","J")
df[,cols] <- 1
df

# So I can just create a vector of all the plate+well combos 
# and initialize them to 0 simultaneously

output <- data.frame(original_seq=seqtab$X)
output[,c("f_barcode","r_barcode","well","trimmed_seq")] <- ""

plate_rows <- c("A","B","C","D","E","F","G","H")
plate_cols <- 1:12
wells <- as.vector(t(outer(plate_rows, plate_cols, paste, sep="")))
plate_names <- colnames(seqtab)[-1]
plate_well_combos <- as.vector(t(outer(plate_names, wells, paste, sep="_")))

output[,plate_well_combos] <- 0

# Main loop
for (i in 1:nrow(seqtab)) {
  seq <- seqtab$X[i]
  j <- which(startsWith(seq, barcodes$F) & endsWith(seq, barcodes$R))
  
  if (length(j) == 0) {
    output[i,"f_barcode"] <- NA
    output[i,"r_barcode"] <- NA
    output[i,"well"] <- NA
    next
  } 
  
  f_barcode <- barcodes[j,1]
  r_barcode <- barcodes[j,2]
  well <- barcodes[j,3]
  
  output[i,"f_barcode"] <- f_barcode
  output[i,"r_barcode"] <- r_barcode
  output[i,"well"] <- well
  
  # Trim sequence
  len_front <- nchar(f_primer) + nchar(f_barcode)
  len_back <- nchar(r_primer) + nchar(r_barcode)
  len_seq <- nchar(seq)
  
  output[i,"trimmed_seq"] <- substr(seq, len_front+1, len_seq-len_back)
  
  # Populate counts
  # For current row, find which plates have non-zero counts
  # For each plate with a non-zero count, assign total to that plate + this well
}
