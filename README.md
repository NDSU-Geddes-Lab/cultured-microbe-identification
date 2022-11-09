# cultured-microbe-identification
### Identify cultured microbes and calculate their purity

### Prerequisites:
	1. R (4.1.0 and above)
	2. Rstudio
	3. Python (3.9 and above) 
	4. R packages:
		1. dada2 (It is always a good idea to make sure that the dada2 package is already installed and running)
		2. tidyverse
		3. reticulate


### Organization of the data:
	1. Create a new main folder to store your input, intermediate and output files.
	2. In this main folder,  create two folders input and output. Store your fastq files in the input folder and 
       keep the output folder for the final result files.
	3. Store the `BC_to_well2.csv file` in the main folder.
	4. Copy and store the `main_workflow.Rmd` and `process_dada2_seqtab.py` files in the main folder.
	5. Go to the Rstudio, Click on the drop down menu named as "Project: (None)" and select "New Project". 
       Click on the "Existing Directory". Click on the "Browse".Select the main folder created in the beginning.
	6. After doing this successfully, there should be input folder, output folder, `main_workflow.Rmd` and 
       `process_dada2_seqtab.py` in the right bottom panel on the Rstudio.


### Execution of the scripts:
	1. `main_workflow.Rmd` is a R notebook with multiple code chunks.
	2. First code chunk is for loading the required R packages (dada2, tidyverse and reticulate)
	3. The second code chunk is the dada2 workflow part. It will output a couple of things below the code chunk:
		1. List of the fastq files
		2. Quality plots for the fastq files
	4. The second code chunk will also output intermediate output files:
		1. A folder named as "filtered" in the input folder containing filtered fastq files
		2. A file named *"input_seqtab.csv"*
	5. The third code chunk loads pyhton3 for the next step.
	6. The fourth code chunk runs the python script that process each sequence in the *"input_seqtab.csv"* file.
		1. For this part, user will provide appropriate arguments to the parameters:
		2. The forward primer sequence (--fprimer), reverse primer sequence (--rprimer), 
            path to the input_seqtab.cvs (--seqtab-path), barcode csv file (--barcode-path) and 
            processed_data_output.csv (--outpath-path).
		3. Once the all arguments are in place, run the chunk
	7. By this time we have an intermediate output in a csv file called `processed_data_output.csv` which holds following 
       information for each sequence:
		1. Original sequences
		2. Its identified forward barcode
		3. Its identified reverse barcode
		4. Identified well name
		5. Counts of the sequences in every well of every plate (zero count in the well except the identified well is expected)
	8. Next code chunk collapses the processed data such that we get a table with all the unique clean sequences with 
       their counts in all the wells in all the plates (960 columns named as Plate*_well*).
	9. Now we calculate purity percentage of each unique clean sequence. For each unique clean sequences, its counts 
       are sorted from highest to lowest. It divide each count by the total number of counts for that Plate*_well* and 
       multiplies by 100.
 
