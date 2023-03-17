'''
import os 

#step one - download FASTQ files associated with the HCMV transcriptomes
print('step one -- downloading transcriptomes') 


#download HCMV transcriptomes from SRA by first assigning the site names to the associated SRR ID
SRR5660030 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030' 
SRR5660033 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033' 
SRR5660044 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044'
SRR5660045 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045'
 
#then using wget to download from the associated link 
os.system('wget ' + SRR5660030)
os.system('wget ' + SRR5660033)
os.system('wget ' + SRR5660044)
os.system('wget ' + SRR5660045) 

#use fastq-dump to donwload FASTQs from SRAs
os.system('fastq-dump -I --split-files SRR5660030') 
os.system('fastq-dump -I --split-files SRR5660033')
os.system('fastq-dump -I --split-files SRR5660044')
os.system('fastq-dump -I --split-files SRR5660045')

print('finished downloading HCMV transcriptomes... ') 

os.mkdir('SampleData') 
'''
'''
import os

SRRfastqs = ['SRR5660030_1.fastq', 'SRR5660030_2.fastq', 'SRR5660033_1.fastq', 'SRR5660033_2.fastq', 'SRR5660044_1.fastq', 'SRR5660044_2.fastq', 'SRR5660045_1.fastq',  'SRR5660045_2.fastq']

samplesList = []
for i in SRRfastqs: 
    samplesList.append('cat ' + i + ' | head -n 4 > SampleData/' + i ) 
print(samplesList[1])

test = samplesList[1]

os.system(test) 
'''

#Determine if the user wants to run on sample data or full data: 
import os
#if running on sample, please run: 
os.chdir(~\PipelineProject_Anisa_Nasse\SampleData) 
#and continue, otherwise, leave the line above commented out. 


###########################################################################

#step two - which strains are most similar to the patient samples


#Create an index for HCMV to be used for mapping 
import Bio 
from Bio import Entrez 

print('step two -- building transcriptome index for HCMV using kallisto') 

#retrieving HCMV records from NCBI
Entrez.email = 'anasse@luc.edu' 
handle = Entrez.efetch(db = 'nucleotide', id = ['NC_006273.2'], rettype = 'fasta') 
records = handle.read() 
#print record to an outfile for later use 
HCMVrecord = open('HCMVrecord.txt', 'w')  
HCMVrecord.write(records) 
HCMVrecord.close()



'''
#os.system('echo " Donor 1 (HCMV30) had 2259287 read pairs before Bowtie2 filtering and 1480440 read pairs after" >> log ')
#os.system('git status') 
#os.system('git add log') 
#os.system('git commit -m "donor commits"') 
#os.system('git push') 
'''

'''
#Use Bowtie2 to create an index for HCMV (NCBI accession NC_006273.2) 
os.system('bowtie2-build HCMVrecord.txt HCMV')


#building .sam files
# bowtie2 calls the program required
# --quiet prints nothing to the screen except errors
# -x specifies the name of the index
# -1 and -2 are the matching input fastq files
# -S names the output file
print('starting bowtie2 .sam') 
os.system('bowtie2 --quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV30mapped.sam')
print('finished 30') 
os.system('bowtie2 --quiet -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S HCMV33mapped.sam')
print('finished 33') 
os.system('bowtie2 --quiet -x HCMV -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S HCMV44mapped.sam')
print('finished 44') 
os.system('bowtie2 --quiet -x HCMV -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMV45mapped.sam')
print('finished 45')

 
#save only the reads that map to the HCMV index for use in SPAdes assembly
# --al-conc-gz writes paired-end reads that align in the file
print('save only the reads that map to the HCMV index for use in SPAdes assembly') 
os.system('bowtie2 -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV30mapped.sam --al-conc-gz HCMV30_mapped_%.fq.gz')
print('starting 33') 
os.system('bowtie2 -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV30mapped.sam --al-conc-gz HCMV30_mapped_%.fq.gz')
print('starting 33') 
os.system('bowtie2 -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S HCMV33mapped.sam --al-conc-gz HCMV33_mapped_%.fq.gz')
print('starting 44') 
os.system('bowtie2 -x HCMV -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S HCMV44mapped.sam --al-conc-gz HCMV44_mapped_%.fq.gz')
print('starting 45') 
os.system('bowtie2 -x HCMV -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMV45mapped.sam --al-conc-gz HCMV45_mapped_%.fq.gz')


############################################################################


#step 3 - using the Bowtie2 output reads, assemble all four transcriptomes together to produce one assembly via SPAdes. Write SPAdes commany used to the log file.

#spades assembly 

os.system('spades.py -k 77,99,127 -t 2 --only-assembler --pe-1 1 HCMV30_mapped_1.fq.gz --pe-2 1 HCMV30_mapped_2.fq.gz --pe-1 2 HCMV33_mapped_1.fq.gz --pe-2 2 HCMV33_mapped_2.fq.gz --pe-1 3 HCMV44_mapped_1.fq.gz --pe-2 3 HCMV44_mapped_2.fq.gz --pe-1 4 HCMV45_mapped_1.fq.gz --pe-2 4 HCMV45_mapped_2.fq.gz -o HCMV2-SRR_assembly') 


#3. write python code to calculate the number of contigs with a length > 100 and write the number to the log file 

#finding the longest contig 
import Bio
from Bio import SeqIO

seqs = []
with open('HCMV2-SRR_assembly/contigs.fasta') as f:
    for record in SeqIO.parse(f, 'fasta'):
#append each sequence in record to the list seqs 
        seqs.append(record.seq)
        longestContig = (max(seqs, key = len))
#the longest contig is the longest sequence in seqs, name it longestContig
        #print(longestContig) 
with open('longestContig.txt', 'w') as file: 
#write the longestContig to the outfile and save for later use 
    file.write(str(longestContig)) 
file.close()


#write code to calculate the length of the assembly and write to the log file
import Bio
from Bio import SeqIO

j = []
with open('HCMV2-SRR_assembly/contigs.fasta') as f:
    for record in SeqIO.parse(f, 'fasta'):
        j.append(len(record.seq))
#append the length of each sequence in record to the list j 
        count = 0
        for i in j:
            if i > 1000:
#for each sequence in j, if the length is below 1000, add one to count  
                count += 1
                greatContigs = []
            greatContigs.append(i)
#add the counts to the list 

with open('contigCounts.txt', 'w') as outfile:
    n = (sum(greatContigs))
#count the number of contigs and write to out file 
    outfile.write('There are ' + str(n) + ' bp in the assembly' + '\n')
    counts  = (str(count))
#print the counts to the outfile 
    outfile.write('There are ' + counts + ' > 1000 bp in the assembly')
    print('bp counts: ', sum(greatContigs))
    print('contig counts: ', count)


############################################################################


#step 4. Does your assembly align with other virus strains? 


#bringing in the nucleotide records associated with Betaherpesvirinae from NCBI
#esearch: interprets querys from command line 
# -db: specifies the database searched 
# -query: tells esearch what you're looking for 
#efetch: outputs esearch in a specified form 
# -format: tells the program the file format to output 
# >: output file name 
os.system('esearch -db Nucleotide -query ""Betaherpesvirinae"[Organism] OR Betaherpesvirinae[All Fields]" | efetch -format fasta > BetaherpesvirinaeRecords.txt')


#making a local database using blast 
#makeblastdb: calls the program requuired 
# -in: name of infile 
# -out: name of outfile 
# -title: name of the local database 
# -dbtype: type of database used 
os.system('makeblastdb -in BetaherpesvirinaeRecords.txt  -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl')  

import os 

#step one - download FASTQ files associated with the HCMV transcriptomes
print('step one -- downloading transcriptomes') 


#download HCMV transcriptomes from SRA by first assigning the site names to the associated SRR ID
SRR5660030 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030' 
SRR5660033 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033' 
SRR5660044 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044'
SRR5660045 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045'
 
#then using wget to download from the associated link 
os.system('wget ' + SRR5660030)

'''
