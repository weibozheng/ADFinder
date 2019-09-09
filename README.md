# ADFinder
ADFinder is a software package that identifies programmed DNA eliminations (PDEs) using large-scale mapping of unlimited NGS datasets. ADFinder is compatible with NGS genomic reads as the input data and can quantitatively detect all PDEs. By using the “local” alignment mode of bowtie2, ADFinder can directly analyze the clipped reads and record the number of clipped and non-clipped reads at each splice site, then it obtains the deletion segments at each splice site and calculate the possibility of each splice (deletion) event. 

# Usage
This software is written in Python3.  Bowtie2 software is required and should be available in PATH.
**Required module:** Biopython, threading, queue.  
You can install these modules by: "pip install bio (or threading/queue)".  
**Usage:** "ADFinder_singleT.py -p <number of threads> -x <genome reference in fasta> [[--s1 <paired-end reads 1_L,2_L> --s2 <paired-end reads 1_R,2_R>]&[-U <Unpaired reads 1,2,3>]]|[-S sam file] -o <output folder>"  
**One example: **  
  "python ADFinder_singleT.py -p 4 -x f1.fasta --s1 f1.left.fq.gz,f2.left.fq.gz,... --s2 f1.right.fq.gz,f2.left.fq.gz,... -U unpaired.fq.gz,... -o folder_output  
or:  
"python ADFinder_singleT.py -p 4 -x f1.fasta -S mapping.sam -o folder_output"  

**output files**  
"output.deletion.tab": deletion segments;  
"output.insertion.tab": insertion segments;  
"output.splice": detailed records of all splices;  
"output.uncompleted.splice": uncompleted splices, containing insertions longer than 300 bp.  
