# cultured-microbe-identification
### Identify cultured microbes and calculate their purity


### To download the scripts from GitHub:

- Go to the cultured-microbiome-identification repository
- Click on the **Code** 
- From the dropdown list, select **Download Zip** option
- Move this zipped repository to the repository that holds your input fastq files


### Prerequisites:

	1. R (4.1.0 and above)
	2. Rstudio
	3. Python (3.9 and above) 
	4. R packages:
		1. dada2 (It is always a good idea to make sure that the dada2 package is already installed and running)
		2. tidyverse
		3. reticulate


##### Preparation for python3 environment on Mac:

- Open terminal on Mac
- Change directories to the directory where you downloaded the scripts
- Run `pip3 install virtualenv` (if virtualenv not installed)
- Next, run `python3 -m venv python_env` (name your python virtual environment, here we named it **python_env**) 
- Lastly, run `source python_env/bin/activate` to activate the environment


##### Preparation for python3 environment on Windows:



### Organization of the data:

- Create a new main folder to store your input, intermediate and output files. 
- In this main folder,  create two folders input and output. Store your fastq files in the input folder and 
keep the output folder for the final result files.
- Store the `BC_to_well2.csv` file in the main folder.
- Copy and store the `main_workflow.Rmd` and `process_dada2_seqtab.py` files in the main folder.
- Go to the Rstudio, Click on the drop down menu named as "Project: (None)" and select "New Project".
- Click on the "Existing Directory". Click on the "Browse". Select the main folder created in the beginning.
- After doing this successfully, there should be input folder, output folder, `main_workflow.Rmd` and 
`process_dada2_seqtab.py` in the right bottom panel on the Rstudio.


### Execution of the scripts:

- `main_workflow.Rmd` is a R notebook with multiple code chunks.
- First code chunk is for loading the required R packages (dada2, tidyverse and reticulate)
- The second code chunk is the dada2 workflow part. It will output a couple of things below the code chunk:
  1. List of the fastq files
  2. Quality plots for the fastq files
- The second code chunk will also output intermediate output files:
  1. A folder named as "filtered" in the input folder containing filtered fastq files
  2. A file named **"input_seqtab.csv"**
- The third code chunk loads python3 for the next step.
- The fourth code chunk runs the python script that process each sequence in the **"input_seqtab.csv"** file. 
  - For this part, user will provide appropriate arguments to the parameters:
    1. The forward primer sequence (--fprimer), reverse primer sequence (--rprimer), 
    path to the input_seqtab.cvs (--seqtab-path), barcode csv file (--barcode-path) and 
    processed_data_output.csv (--outpath-path).
    2. Once the all arguments are in place, run the chunk
- By this time we have an intermediate output in a csv file called `processed_data_output.csv` which holds following 
information for each sequence:
          1. Original sequences
          2. Its identified forward barcode
          3. Its identified reverse barcode
          4. Identified well name
          5. Counts of the sequences in every well of every plate (zero count in the well except the identified well is expected)
- Next code chunk collapses the processed data such that we get a table with all the unique clean sequences with their 
counts in all the wells in all the plates (960 columns named as Plate*_well*).
- Now we calculate purity percentage of each unique clean sequence. For each unique clean sequences, its counts are 
sorted from highest to lowest. It divides each count by the total number of counts for that Plate*_well* and multiplies 
by 100.
 
