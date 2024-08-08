# Cultured Microbe ID workflow

Identify cultured microbes and calculate their purity

## Workflow overview

The following is a high-level overview of the Cultured Microbe ID workflow.

1. Raw reads are processed using a basic [DADA2](https://benjjneb.github.io/dada2/) processing workflow resulting in a `seqtab` object, where rownames correspond to unique sequences and column names correspond to culture plate IDs. Each entry in the `seqtab` is a count representing the number of times a particular sequence appeared in a particular plate.
2. For each sequence in the `seqtab`,the pair of forward and reverse barcodes from the `BC_to_well2.csv` plate map are identified, which determines the specific plate well the sequence came from. Barcodes and primers are then trimmed and counts are recorded in that plate+well column.
3. From the set of trimmed sequences in step 2, unique trimmed sequences are identified and counts for each unique sequence are summed across all rows, resulting in a data frame containing count information for each sequence across all plates and wells.
4. The top *n* plate wells for each trimmed sequence are identified by sorting the columns for that row. For the top wells, a percetage purity is calculated by dividing the count by the total of all counts for that plate well and multiplying by 100.
5. Lastly, taxonomy is assigned to each unique sequence, and a results are written to a CSV file. Now for each unique sequence, we have the following information:
  - ASV
  - Taxonomy
  - Top *n* counts of the ASV and the purity of the ASV in the well

## Prerequisites and dependencies

In order to run this workflow, you will need to install the following:

	1. [R language](https://www.r-project.org/) (4.1.0 and above)
	2. [RStudio IDE](https://posit.co/download/rstudio-desktop/) (if you intend to run the interactive RMarkdown notebook)
	3. Add-on R libraries

```r
install.packages(c("tidyverse","argparser","BiocManager"))
BiocManager::install("dada2")
```

## Downloading the code

You can download the source code for this workflow by cloning the repository like so:

```bash
git clone https://github.com/NDSU-Geddes-Lab/cultured-microbe-identification.git
```

Or you can download this repository manually by following the steps below.

- Click on the **Code**.
- From the dropdown list, select **Download Zip** option.
- Move this zipped repository to the repository that holds your input fastq files and unzip it.

## Preparing the analysis

- After downloading and extracting the source code from above step, you should have a directory called `cultured-microbe-identification`.
- In the `cultured-microbe-identification` folder, create two sub-folders `input` and `output`.
- Place your fastq files to be processed in the `input` sub-folder.

### Execution

Now that you have all the pre-requisites setup and you understand high-level workflow, please open `main_workflow.Rmd` file in Rstudio.

`main_workflow.Rmd` file contains all step-by-step instructions on how to run each code chunk or section and if you are required to modify any value in any chunk.
