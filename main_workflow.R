library(tidyverse)
library(dada2)

# Path to Silva database
silva_db_path <- "./db/silva_nr99_v138.1_train_set.fa.gz"

# Path to the folder where all the input fastq files are stored
path <- "./input"

# List all the files present in `path` variable
list.files(path)

#Sort files to ensure forward/reverse are in the same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

#extract sample names, assuming filenames have format:
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Now we visualize the quality profile of the reverse reads:
plotQualityProfile(fnFs[3:6])
plotQualityProfile(fnRs[3:6])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Get the total number of input fastq files
length(fnFs)
length(fnRs)

# Check if there are any duplicate input fastq files
any(duplicated(c(fnFs, fnRs)))
any(duplicated(c(filtFs, filtRs)))

# Check of total number of filtered fastq files
length(filtFs)
length(filtRs)

#Trim based on quality plots, in this example we cut for the forward reads to 240bp and the reverse reads to 200bp
out <- filterAndTrim(fnFs, 
                    filtFs, 
                    fnRs, 
                    filtRs, 
                    maxN=0, 
                    truncLen=c(240,200),
                    maxEE=c(2,2), 
                    truncQ=2, 
                    rm.phix=TRUE, 
                    compress=TRUE, 
                    multithread=TRUE)  # On Windows set multithread=FALSE
head(out)

#Plot again to see if the trimming worked
plotQualityProfile(filtFs[3:6])
plotQualityProfile(filtRs[3:6])

table(file.exists(filtFs))
table(file.exists(filtRs))

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plot error model
plotErrors(errF, nominalQ=TRUE)

#Apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

#Merge paires reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[3]])
head(mergers[[10]])

#Construct sequence table to ASV
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
# the amplicon should be around 299 to 303 bp long
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 299:303]
table(nchar(getSequences(seqtab2)))

# Transpose seqtab for processing barcodes
seqtab2 <- t(seqtab2)

# Load barcodes from plate map
barcodes <- read.csv("BC_to_well2.csv")

# Specify forward and reverse primers
f_primer <- "GTGCCAGCMGCCGCGGTAA"
r_primer <- "GACTACHVGGGTATCTAATCC"

# Create output data frame for barcode matches
processed_data <- data.frame(original_seq=rownames(seqtab2))
processed_data[,c("f_barcode","r_barcode","well","trimmed_seq")] <- ""

# Create list of plate+well combos
plate_rows <- c("A","B","C","D","E","F","G","H")
plate_cols <- 1:12
wells <- as.vector(t(outer(plate_rows, plate_cols, paste, sep="")))
plates <- colnames(seqtab2)
plate_well_combos <- as.vector(t(outer(plates, wells, paste, sep="_")))

# Add plate+well combos to data frame and initialize to 0
processed_data[,plate_well_combos] <- 0

# For each sequence in the seqtab, determine the plate well
# by identifying the matching forward and reverse barcodes.
# Then assign the count for that sequence to the appropriate plate well.
for (i in 1:nrow(seqtab2)) {
  seq <- rownames(seqtab2)[i]
  
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
  processed_data[i,"f_barcode"] <- f_barcode
  processed_data[i,"r_barcode"] <- r_barcode
  processed_data[i,"well"] <- well
  
  # Trim sequence
  # (regardless of whether there is a match or not, since that's what Komal did)
  len_front <- nchar(f_primer) + nchar(f_barcode)
  len_back <- nchar(r_primer) + nchar(r_barcode)
  len_seq <- nchar(seq)
  
  # Add trimmed sequence to output
  processed_data[i,"trimmed_seq"] <- substr(seq, len_front+1, len_seq-len_back)
  
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
    seq_count <- seqtab2[i,plate]
    if (seq_count > 0) {
      plate_well <- paste(plate, well, sep="_")
      processed_data[i,plate_well] <- seq_count
    }
  }
}

# Extract trimmed seqs and their count information
trimmed_seqs_data <- processed_data %>% 
  select("trimmed_seq", all_of(plate_well_combos))

# Collapse (sum) plate-well values for each unique trimmed 
unique_trimmed_seqs <- aggregate(. ~ trimmed_seq, data=trimmed_seqs_data, FUN=sum)

rownames(unique_trimmed_seqs) <- unique_trimmed_seqs$trimmed_seq
unique_trimmed_seqs$trimmed_seq <- NULL

# Additional functions needed for purity calculations
source("filter_purity.R")

# Identify wells with top n counts
final_df <- process_each_sequence(unique_trimmed_seqs, 5)

# Taxonomy
tax <- get_taxonomy(final_df$ASV, silva_db_path)

# Combine taxonomy with purity information
fd <- final_df
rownames(fd) <- fd$ASV
fd_with_taxa <- merge(tax, fd, by = 0, all = TRUE)
names(fd_with_taxa)[names(fd_with_taxa) == "Row.names"] <- "ASV"

# Write final data frame into a csv file
write.csv(fd_with_taxa, file = "./output/asv_analysis_results.csv")

