# Lab Notebook

**Use terminal**

**To login: ssh lcaicedo@pbio381.uvm.edu**

**Password: (mine)**

01/31/2017

Definitions:

$ means waiting for me to issue a command

root is the "God" of the machine

<u>Basic commands:</u>

ls = list all files

top = shows all users and what they are doing

ll = "list in long format" show current contents of folders (shows size)

q = get out of top

~/ = go to my home directory

pwd = tells me where I am, e.g. /users/l|c/lcaicedo (l|c means name broken into subdirectories)

mkdir = make a new folder, e.g.. mkdir mydata

**cd foldername** = move within a directory

cd /data/ = goes to most basal shared space (if i need to make a copy for myself)

cd project-data/ = look for directory nested within the one we are looking for.

cp = copy



mv = to move but also can be used to rename files 

​	mv samp_fullfilename newname

screen = opens a new space (window) a d type command (eg. bash bwaaln.sh) then cntl + a + d (keys) —— this will detach

- to come back: screen-r (re-attach)

02/06/2017

cd /data/project_data/fastq

my sample. 07 _ 5-11_ _S_ _4 _R1.fq.gz: 07 = individual, 5-11 = date, S = Sick, 1= stage1 (sampled on R1= right read bc they are paired reads) + R2 (left read)

zcat - program to open zipp files

- zcat FILENAME | head



Quality (Phred) scores: probablity of being wrong

- asci format: double digit value of quality score be represented by one letter or symbol (high prob)
- ig Q >= 30, equivalent to 0.1% chance of error



FastQC= program to look at quality scores systematically, clean the data file

- first: cp trim_example.sh ~/scripts/
- cd script/
- vim trim_example.sh~/scripts/
- i = to edit

Trimming: Trimmomatic program

/data/project_data/fastq/yoursample_R1.fq.gz \

/data/project_data/fastq/yoursample_R2.fq.gz \

/data/project_data/fastq/cleanreads/cleanreads/"samp_R1_clean_paired.fq" \

/data/project_data/fastq/cleanreads/cleanreads/"samp_R1_clean_unpaired.fq" \

/data/project_data/fastq/cleanreads/cleanreads/"samp_R2_clean_paired.fq" \

/data/project_data/fastq/cleanreads/cleanreads/"samp_R2_clean_unpaired.fq" \

- esc key
- :w!cd = to write (save changes)
- :q! = to quit
- scp = direct html file to you computer 
- cnt + c (keys) = kill a process

________

<u>02/08/2017: Transcriptomics 2</u>

Goals:

- finish cleaning (trim)
- fastqc (vis)
- make table of number of reads
- design assembly tests
- start assemblies
- evaluate assembly

Learning Objectives (skills)

- scripts
  - paths
    - program
    - input
    - output
  - filenames
    - in
    - out
    - moving directories, files
  - scp
    - move from server to pc
  - Executing scripts

1) Finish cleaning: vim trim_example

2) move to my directory: 

cd .. = to go to previous directory

.trim_example.sh or bash trim_example.sh

fastqc 07_5-11_S_4_R*.fq

on another tab in terminal:

scp lcaicedo@pbio381.uvm.edu:data/project_data/fastq/07_5-11_S_R*.html .

now I have html files on my computer to see trim 

What affects my assembly: the lenght of the read

- use of paired reads

________

02/13/2017

A way of cleaning up is to predict open reading frames: ORF= ___ to say we want complete reads 

QC: rudimentary quantitative control, using N50… count of number of reads

BUSCO = to see how many single copy orthologs are present "provides quantitative measures for the assessment of genome assembly, gene set, and transcriptome completeness, based on evolutionarily-informed expectations of gene content from near-universal single-copy orthologs selected from" - bunco.ezlab.org



*We have the reference transcpritome (1st part of tutorial explains how it was produced)* 



we take individual sample and map to the assembly

07 _ 5-11_ _S_ _4 _R1.fq.gz

variables: myLeft, myRight, myShort

echo= to test that I correctly designed this. 

bwa (program)

​	bwa aln : gives me glossary of symbols. 

to run: bash bwaaln_lc.sh

to run and copy output a text (log file) into scripts folder

when done I should have a .sam (sequence alignment file)

​	FLAG = how well it mapped 

---

02/15/2017

~/ = to go to my directory

cd scripts

Shows what is in the folder: 

total 1932280

-rw-r--r--. 1 lcaicedo users 1926280714 Feb 13 11:14 07_5-11_S_4_bwaaln.sam

-rw-r--r--. 1 lcaicedo users   25887148 Feb 13 11:11 07_5-11_S_4_R1.fq.gz_left_clean_paired.fq.sai

-rw-r--r--. 1 lcaicedo users   26467948 Feb 13 11:13 07_5-11_S_4_R2.fq.gz_right_clean_paired.fq.sai

-rwxr-xr-x. 1 lcaicedo users        903 Feb 13 11:07 bwaaln_lc.sh

-rw-r--r--. 1 lcaicedo users        835 Feb  8 10:44 cd