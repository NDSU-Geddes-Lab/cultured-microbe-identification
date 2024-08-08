# Cultured Microbe ID workflow

Identify cultured microbes and calculate their purity

## Workflow overview

The following is a high-level overview of the Cultured Microbe ID workflow.

1. Raw reads are processed using a basic [DADA2](https://benjjneb.github.io/dada2/) processing workflow resulting in a `seqtab` object, where rownames correspond to unique sequences and column names correspond to culture plate IDs. Each entry in the `seqtab` is a count representing the number of times a particular sequence appeared in a particular plate.
2. For each sequence in the `seqtab`,the pair of forward and reverse barcodes from the `BC_to_well2.csv` plate map are identified, which determines the specific plate well the sequence came from. Barcodes and primers are then trimmed and counts are recorded in that plate+well column.
3. From the set of trimmed sequences above, unique trimmed sequences are identified and counts for each unique sequence are summed across all rows, resulting in a data frame containing count information for each sequence across all plates and wells.
4. The top *n* plate wells for each trimmed sequence are identified by sorting the columns for that row. For the top wells, a percetage purity is calculated by dividing the count by the total of all counts for that plate well and multiplying by 100.
5. Lastly, taxonomy is assigned to each unique sequence, and a results are written to a CSV file. Now for each unique sequence, we have the following information:
  - ASV
  - Taxonomy
  - Plates and wells containing the top *n* counts for each unique ASV and the purity of the ASV in the well

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
- In the `cultured-microbe-identification` folder, create three sub-folders: `input`, `output`, and `db`.
- Place your fastq files to be processed in the `input` sub-folder.
- Place your Silva taxonomy database in the `db` folder.

*Note: These are just the defaults. Paths to input, output, and database folders can be changed manually in the RMarkdown notebook or via command line arguments in the CLI script.*

## Running the workflow

There are two ways to run this workflow:

1. Executing code chunks interactively in RStudio using the `main_workflow.Rmd` notebook.
2. Executing the entire workflow via the command line using the `main_workflow.R` CLI script.

Each method is described in detail below.

### 1. Interactive execution via RStudio and RMarkdown

To run the workflow interactively:

1. Use the "Files" menu in the bottom-right pane of RStudio to navigate to the `cultured-microbe-identification` folder you prepared above.
2. Open the `main_workflow.Rmd` file in the RStudio editor.
3. Execute each code chunk, in order, by clicking the green triangle (i.e., "play" button) in the upper-right corner of the chunk.
4. Outputs will be available in the `ouptut` subdirectory you created above. Plots and other diagnostic messages will appear directly in the notebook.

### 2. Command line execution using the CLI script

To run the workflow from the command line:

1. Open a terminal and navigate to the `cultured-microbe-identification` folder created above.
2. To see a list of options for the CLI script, you can run it with the `-h` flag, like so:

```bash
Rscript main_workflow.R -h
```

You should then see a help menu like this:

```
usage: main_workflow.R [--] [--help] [--db DB] [--barcodes BARCODES]
       [--fwd FWD] [--rev REV] [--hits HITS] [--outdir OUTDIR]
       fastq_dir

Cultured Microbe ID

positional arguments:
  fastq_dir       Directory containing input sequences (.fastq.gz)

flags:
  -h, --help      show this help message and exit

optional arguments:
  -d, --db        Path to taxonomy database [default:
                  ./db/silva_nr99_v138.1_train_set.fa.gz]
  -b, --barcodes  Path to barcode plate map (.csv) [default:
                  ./BC_to_well2.csv]
  -f, --fwd       Forward primer [default: GTGCCAGCMGCCGCGGTAA]
  -r, --rev       Reverse primer [default: GACTACHVGGGTATCTAATCC]
  --hits          Number of hits to report (top n wells) [default: 5]
  -o, --outdir    Output directory [default: ./output]

```

3. To run the workflow, run the script the same way, but provide the location to the input directory containing your sequence reads.

```bash
Rscript main_workflow.R input
```

(Output omitted)

As indicated in the help menu above, the taxonomy database, barcode plate well map, primers, number of hits to report, and output directory can all be changed if desired.
