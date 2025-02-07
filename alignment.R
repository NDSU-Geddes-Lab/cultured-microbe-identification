# Load necessary libraries
library(DECIPHER)
library(phangorn)
library(ape)

checkValidPhylo(treeNJ)

#Change depending on where your files are
path <- ""
list.files(path)
setwd(path)

# Read sequences from a FASTA file
fasta_file <- "corn_HTPC_Slur.fasta" # Replace with the path to your FASTA file
sequences <- readDNAStringSet(fasta_file)

# Run Sequence Alignment (MSA) using DECIPHER
alignment <- AlignSeqs(sequences, anchor = NA)

# Convert sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")

# Create a distance matrix
dm <- dist.ml(phang.align)

# Perform Neighbor Joining to construct an initial tree
treeNJ <- NJ(dm)
treeNJ <- root(treeNJ, outgroup = "C1_Slur", resolve.root = TRUE)

# Fit a maximum likelihood model on the tree
fit <- pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

# Extract the final optimized tree
final_tree <- fitGTR$tree

# Export the tree to a Newick file
write.tree(final_tree, file = "final_tree.newick")

# Optionally, plot the tree to visualize
plot(final_tree, main = "Phylogenetic Tree from FASTA Alignment")
