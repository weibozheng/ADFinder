# ADFinder
Searching for accurate deletion/insertion segments using NGS data

# Usage
This software is written in Python3.  
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
