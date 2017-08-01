# allele-wrangler
For Wrangling HLA Alleles

Python 2.7

Reads are self-aligned to create a consensus sequence.

This requires Anaconda to run.

## To configure Anaconda
Anaconda uses separate environments to run your programs in.  
Install Anaconda for python 2.7.  
https://www.continuum.io/downloads  
To set up the environment in anaconda:  

Linux/Mac:  
```
conda create --name minionvironment biopython six pycurl pysam
source activate minionvironment  
pip install pyinstaller packaging matplotlib
source deactivate  
```  
Windows:  
```  
conda create --name minionvironment biopython six pycurl pywin32 pysam  
call activate minionvironment && pip install pyinstaller packaging matplotlib && call deactivate  
```

You must also have clustalo installed for this program to generate initial consensus.
'clustalo'
i used apt-get to install that on ubuntu.  On Windows:
http://www.clustal.org/omega/

It uses bw aligner for iterating the reads.
Also used ubuntu to install that.  'bwa'


Future: can we split groups of reads for heterozygous alleles, and then assemble?

Future: Can we HLA allele call?


