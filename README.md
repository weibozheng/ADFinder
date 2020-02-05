# ADFinder

Programmed DNA elimination plays a crucial role in the transitions between germline and somatic genomes in diverse organisms ranging from unicellular ciliates to multicellular nematodes. However, software specific to the detection of DNA splicing events is scarce. Here we develop ADFinder, an efficient detector of programmed DNA eliminations using NGS high-throughput sequencing data. ADFinder can predict programmed DNA eliminations with relatively low sequencing coverage, detect multiple alternative splicing forms in the same genomic location, and calculate the frequency for each splicing event. This software would facilitate research of programmed DNA eliminations and all down-stream analyses.

# Description
This software is written in Python3.  Bowtie2 software is required and should be available in PATH.
#Required python3 module
Biopython.  
You can install these modules by: "pip install bio ".  

## Usage
ADFinder.py **-x** \<genome reference in fasta format\> [**--s1** \<paired-end reads 1\> **--s2** \<paired-end reads 2\>]
|[**-U** \<unpaired reads\>]\|\[**-S** \<sam file\>] **-o** \<output folder\>

**One example:**  
  "python ADFinder_singleT.py -p 4 -x f1.fasta --s1 f1.left.fq.gz,f2.left.fq.gz,... --s2 f1.right.fq.gz,f2.left.fq.gz,... -U unpaired.fq.gz,... -o folder_output  
or:  
"python ADFinder_singleT.py -p 4 -x f1.fasta -S mapping.sam -o folder_output"  

## output files 
"output.deletion.tab": a nine-column file recording all deletion events;

### format of output file

|tag|chromosome|start|end|sequence|length|splicing_depth|non_splicing_depth|frequency|
|---|---|---|---|---|---|---|---|---|
|D1|OXYTRI_MIC_78915|6767|6785|GGATAATATATTTTTATAT|19|15|1.500000|0.909091|
|D2|OXYTRI_MIC_78915|5785|5809|AATCCTTACTTTAAAATTTTTTATT|25|22|3.100000|0.876494|
|D3|OXYTRI_MIC_78915|4476|4505|AAATTTAATTAAATTTATATGTTTAATTTT|30|7|1.800000|0.795455|


