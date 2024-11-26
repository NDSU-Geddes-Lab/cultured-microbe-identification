suppressPackageStartupMessages({
  library(tidyverse)
  library(dada2)
  library(argparser)
  library(logger)
})

# Process command line args
parser <- arg_parser("Cultured Microbe ID", hide.opts = TRUE)

parser <- add_argument(parser, "fastq_dir", 
                       help = "Directory containing input sequences (.fastq.gz)")

parser <- add_argument(parser, "--db", 
                       help = "Path to taxonomy database", 
                       default = "./db/silva_nr99_v138.1_train_set.fa.gz")

parser <- add_argument(parser, "--barcodes", 
                       help = "Path to barcode plate map (.csv)", 
                       default = "./BC_to_well2.csv")

parser <- add_argument(parser, "--fwd", 
                       help = "Forward primer", 
                       default = "GTGCCAGCMGCCGCGGTAA")

parser <- add_argument(parser, "--rev", 
                       help = "Reverse primer", 
                       default = "GACTACHVGGGTATCTAATCC")

parser <- add_argument(parser, "--hits", 
                       help = "Number of hits to report (top n wells)", 
                       default = "5")

parser <- add_argument(parser, "--outdir", 
                       help = "Output directory", 
                       default = "./output")

argv <- parse_args(parser)

# Prepare output directory and place to save plots
dir.create(file.path(argv$outdir), showWarnings=FALSE)
pdf(file.path(argv$outdir, "Rplots.pdf"))

# Path to Silva database
silva_db_path <- argv$db

# Path to the folder where all the input fastq files are stored
path <- argv$fastq_dir

#Sort files to ensure forward/reverse are in the same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Check number of forward and reverse reads
if (length(fnFs) != length(fnRs)) {
  log_error("Different numbers of forward and reverse read files. Exiting...")
  quit()
}

#extract sample names, assuming filenames have format:
sample.names.fwd <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names.rev <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

# Make sure F and R reads match up or fail with error.
if (! all(sample.names.fwd == sample.names.rev)) {
  log_error("Forward and reverse sample names do not match. Exiting...")
  quit()
}

sample.names <- sample.names.fwd

# Lastly, check for duplicate sample names
if (any(duplicated(sample.names))) {
  duplicates <- unique(sample.names[duplicated(sample.names)])
  log_error("Duplicate sample names: {duplicates}")
  quit()
}

# Print sample names for user.
log_info("Identified {length(sample.names)} samples: {paste(sample.names, collapse=' ')}")

#Now we visualize the quality profile of the reverse reads:
#plotQualityProfile(fnFs)
#plotQualityProfile(fnRs)

# Place filtered files in filtered/ subdirectory
dir.create(file.path(argv$outdir, "filtered"), showWarnings=FALSE)
filtFs <- file.path(argv$outdir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(argv$outdir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

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

#Plot again to see if the trimming worked
#plotQualityProfile(filtFs)
#plotQualityProfile(filtRs)

# Check to see if any samples were dropped after filtering
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

if (length(filtFs) != length(sample.names)) {
  dropped <- sample.names[!exists]
  log_warn("Dropped {length(dropped)} samples that failed quality filtering: paste(dropped, collapse=' ')")
}

log_info("{length(filtFs)} samples remaining after filtering.")

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plot error model
#plotErrors(errF, nominalQ=TRUE)

#Apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object:
#dadaFs[[1]]
#dadaRs[[1]]

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
barcodes <- read.csv(argv$barcodes)

# Specify forward and reverse primers
f_primer <- argv$fwd
r_primer <- argv$rev

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
final_df <- process_each_sequence(unique_trimmed_seqs, as.integer(argv$hits))

# Taxonomy
tax <- assignTaxonomy(final_df$ASV, silva_db_path, multithread=TRUE)

# Combine taxonomy with purity information
rownames(final_df) <- final_df$ASV
fd_with_taxa <- merge(tax, final_df, by.x = 0, by.y="ASV", all = TRUE)
names(fd_with_taxa)[names(fd_with_taxa) == "Row.names"] <- "ASV"

# Filter out Eukaryotes (if any)
fd_with_taxa <- fd_with_taxa %>% filter(Kingdom != "Eukaryota")

# Write final data frame into a csv file
write.csv(fd_with_taxa, 
          file = file.path(argv$outdir, "asv_analysis_results.csv"))

