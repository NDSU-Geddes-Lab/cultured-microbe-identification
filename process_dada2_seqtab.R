#!/usr/bin/env Rscript

# This script is an R port of Komal Marathe's process_dada2_seqtab.py.
# It does the following:
# 
# 1. Creates an output data frame with the following columns:
#      - Original sequence from the seqtab (original_seq)
#      - Forward and reverse matched barcodes (f_barcode, r_barcode)
#      - Well corresponding to barcode pair (well)
#      - Trimmed sequence with n bases removed, where n = barcode+primer length (trimmed_seq)
#      - A column for each plate+well combination, e.g., "plate1_A1", "plate1_A2", ...
# 2. Then, for each sequence in the seqtab, we do the following:
#      - Identify the forward and reverse barcodes
#      - Trim the primers
#      - Find the corresponding well number in the barcode table
#      - Report the result in the column corresponding to the plate+well combo

seqtab <- read.csv("output/input_seqtab.csv")
# Output from Komal's python code, for comparison
output_exp <- read.csv("output/processed_data_output.csv", row.names=1)
barcodes <- read.csv("BC_to_well2.csv")

f_primer <- "GTGCCAGCMGCCGCGGTAA"
r_primer <- "GACTACHVGGGTATCTAATCC"

output <- data.frame(original_seq=seqtab$X)
output[,c("f_barcode","r_barcode","well","trimmed_seq")] <- ""

plate_rows <- c("A","B","C","D","E","F","G","H")
plate_cols <- 1:12
wells <- as.vector(t(outer(plate_rows, plate_cols, paste, sep="")))
plates <- colnames(seqtab)[-1] # May have to fix index later...
plate_well_combos <- as.vector(t(outer(plates, wells, paste, sep="_")))

output[,plate_well_combos] <- 0

# Main loop
for (i in 1:nrow(seqtab)) {
  seq <- seqtab$X[i]
  
  f_barcode <- ""
  r_barcode <- ""
  well <- "no_well" # To match what Komal has
  
  # Find line in barcode map with matching fwd and rev barcodes
  j <- which(startsWith(seq, barcodes$F) & endsWith(seq, barcodes$R))
  
  # There should be either a single match (length(j)==1) 
  # or no matches (length(j)==0), so if length(j) > 0, set
  # barcodes and well accordingly.
  if (length(j) > 0) {
    f_barcode <- barcodes[j,1]
    r_barcode <- barcodes[j,2]
    well <- barcodes[j,3]
  } 
  
  # Populate corresponding output fields
  output[i,"f_barcode"] <- f_barcode
  output[i,"r_barcode"] <- r_barcode
  output[i,"well"] <- well
  
  # Trim sequence
  # (regardless of whether there is a match or not, since that's what Komal did)
  len_front <- nchar(f_primer) + nchar(f_barcode)
  len_back <- nchar(r_primer) + nchar(r_barcode)
  len_seq <- nchar(seq)
  
  # Add trimmed sequence to output
  output[i,"trimmed_seq"] <- substr(seq, len_front+1, len_seq-len_back)
  
  # If we don't have a barcode match, no point in tabulating plate counts
  # since we don't know what well they come from.
  if (length(j) == 0) {
    next
  }
  
  # Populate counts
  # For current row, find which plates have non-zero counts
  # For each plate with a non-zero count, assign total to
  # the appropriate plate+well combo
  for (plate in plates) {
    seq_count <- seqtab[i,plate]
    if (seq_count > 0) {
      plate_well <- paste(plate, well, sep="_")
      output[i,plate_well] <- seq_count
    }
  }
}

# Validate by comparing R-generated output with output from Komal's python code

# Row names for output_exp start at 0. Fix before comparing.
rownames(output_exp) <- 1:nrow(output_exp)
all.equal(output[,1:1061], output_exp[,1:1061]) #Validation PASSED!
