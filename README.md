# cultured-microbe-identification

## Identify cultured microbes and calculate their purity

Please use this page as the main entrypoint of this workflow for step-by-step instructions on how to use this.

### Prerequisites

Please download and install the following softwares and packages on your laptop/machine.

	1. R (4.1.0 and above)
	2. Rstudio
	3. R packages:
		1. dada2 (It is always a good idea to make sure that the dada2 package is already installed and running)
		2. tidyverse


### Download the scripts from GitHub

Download this repository containing R scripts to do processing in subsequent steps.

- Go to the cultured-microbiome-identification repository
- Click on the **Code** 
- From the dropdown list, select **Download Zip** option
- Move this zipped repository to the repository that holds your input fastq files


### Organization of the data:

- After downloading and extracting the source code from above step, you should have a directory called `cultured-microbe-identification`
- In `cultured-microbe-identification` folder, create two sub-folders `input` and `output`.
- Store your all the fastq files to be processed in the `input` sub-folder and keep the `output` sub-folder for the final result files.
- Store the `BC_to_well2.csv` file in `cultured-microbe-identification` folder.
- Go to the Rstudio, Click on the drop down menu named as "Project: (None)" and select "New Project".
- Click on the "Existing Directory". Click on the "Browse". Select the `cultured-microbe-identification` folder created in the beginning.
- After doing this successfully, there should be `input` sub-folder, `output` sub-folder, `BC_to_well2.csv`, `filter_purity.R`, and `main_workflow.Rmd` in the right bottom panel on the Rstudio.


### Explanation of scripts execution:

This paragraph is only an explanation of what is happening inside the R scripts and is meant for understanding an overview.

Once you read this explanation, then proceed to next section of this page on instructions to actual execution.

- `main_workflow.Rmd` is a R notebook with multiple code chunks.
- First code chunk is for loading the required R packages (dada2, tidyverse)
- The second code chunk is the dada2 workflow part. It will output a couple of things below the code chunk:
  1. List of the fastq files
  2. Quality plots for the fastq files
  3. A folder named as "filtered" in the input folder containing filtered fastq files
- The next part is to map each sequence to a pair of forward and reverse barcodes, trim the barcodes and primers, and attribute sequence counts to the plate+well combo corresponding to the matched barcodes.
- By this time we have an intermediate output, `processed_data_output`, which holds following 
information for each sequence:
  - Original sequences
  - Its identified forward barcode
  - Its identified reverse barcode
  - Identified well name
  - Counts of the sequences in every well of every plate (zero count in the well except the identified well is expected)
- Next code chunk collapses the processed data such that we get a table with all the unique clean sequences with their 
counts in all the wells in all the plates (960 columns named as Plate*_well*).
- Now we calculate purity percentage of each unique clean sequence. For each unique clean sequences, its counts are 
sorted from highest to lowest. It divides each count by the total number of counts for that Plate*_well* and multiplies 
by 100.
 

### Execution

Now that you have all the pre-requisites setup and you understand high-level workflow, please open `main_workflow.Rmd` file in Rstudio.

`main_workflow.Rmd` file contains all step-by-step instructions on how to run each code chunk or section and if you are required to modify any value in any chunk.
