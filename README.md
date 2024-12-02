# Cultured Microbe ID workflow

Identify cultured microbes and calculate their purity.

## Prerequisites and dependencies

In order to run this workflow, you will need to install the following:

 - [R language](https://www.r-project.org/) (4.1.0 and above)
 - Add-on R libraries

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

1. Open a terminal and navigate to the `cultured-microbe-identification` folder created above.
2. To see a list of options for the CLI script, you can run it with the `-h` flag, like so:

```bash
Rscript microbeID.R -h
```

You should then see a help menu like this:

```
usage: microbeID.R [--] [--help] [--db DB] [--barcodes BARCODES]
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
Rscript microbeID.R input
```

(Output omitted)

As indicated in the help menu above, the taxonomy database, barcode plate well map, primers, number of hits to report, and output directory can all be changed if desired.
