# GenomeAssemblyPipeline_Anisa

## Languages and Packages

### Installation for Bash: 

Entrez Direct: E-Utilities on the Unix Command Line for step 4 to bring in the records associated with Betaherpesvirinae from NCBI 

Entrez Direct: https://www.ncbi.nlm.nih.gov/books/NBK179288/ 
```

sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

OR

sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```

Bowtie2: https://www.metagenomics.wiki/tools/bowtie2/install

Bowtie2 GitHub: https://github.com/BenLangmead/bowtie2

```

conda install -c bioconda bowtie2

```

SPAdes: https://cab.spbu.ru/files/release3.12.0/manual.html

SPAdes GitHub: https://github.com/ablab/spades

```
wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz

tar -xzf SPAdes-3.12.0-Linux.tar.gz

cd SPAdes-3.12.0-Linux/bin/
```

### Package Installation for Python: 


import os

import Bio

from Bio import Entrez 

from Bio import SeqIO


## Flags for input files: 

### For Bowtie2 mapping 

bowtie2: calls the program required 

--quiet: prints nothing to the screen except errors 

-x: specifies the name of the index 

-1 and -2: matching input fastq files 

-S: names the output file 


### For saving only the reads that map for use in SPAdes assembly: 

--al-conc-gz: writes paired-end reads that align in the file 


### Running SPAdes assembly: 

spades: calls the program required 

-k: specifies kmer lengths  

-t: number of threads to use

--only-assembler: runs assembly module only

-pe-1 and pe-2: files containing left and right reads for the paired-end libraries 

-o : name of outfile


### Using Esearch: 

esearch: calls the program required 

-db: specifies the database searched 

-query: tell esearch what you're looking for 

efetch: calls the program required, outputs esearch in a specified form 

-format: tells the program the file format to output 



### Making a local database using BLAST: 

makeblastdb: calls the program required 

-in: name of in file 

-out: name of out file 

-title: name the local database 

-dbtype: name of the database used 


## Example Code
To run on the sample code, please download the provided SampleData file onto your linux server, and comment out lines 3 through 26 in the PipelineProject_Anisa_Nasse.py file. 

Or, cd home/anasse/PipelineProject_Anisa_Nasse/SampleData, copy the python file here, comment lines 3 through 26 out and run.  

to run: 

nohup python pipelineproject.py

## Workflow

### Step One: Download FASTQ files associated with the input transcriptomes: 

This step uses wget and NCBI to download required HCMV SRR files for further testing. 

Input Files: 

- SRA downloads from NCBI 

Out Files: 

- SRR files 

- SRR_1 and SRR_2 fastq files 

### Step Two: Finding strains most similar to the patient samples: 

This step uses Entrez to download HCMV records from NCBI, then uses Bowtie2 to create an index for these records. 
The reads that map to each transcriptome are saved for use in step three. 

Input Files: 

- Nucleotide records from NCBI using Entrez

Out Files: 

- Bowtie2 .sam files 

- Bowtie2 mapped reads  

### Step Three: Using output reads from Step Two, assemble the transcriptomes together to produce one assembly: 

SPAdes is used in this step. 

Input File: 

- Mapped reads from step two 

Out File: 

- A directory called 'HCMV2-SRR_assembly' 

### Step Four: Does the assembly from Step Three align with other virus strains? 

Python code is used to find the longest contig produced in Step Four. Blast+ is then used to make a local database for the 
Betaherpesvirinae subfamily, and to use the longest contig to find the top 10 hits most similar to the local database. 

Input File: 

- Betaherpesvirinae records from NCBI

Out Files: 

- longestContig.txt containing the longest contig from the downloaded records file 

- Table containing the top 10 blast command outputs 



