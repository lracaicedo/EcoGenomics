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

cat filename.txt : print screen



______

02/06/2017

```
cd /data/project_data/fastq
```

my sample. 07 _ 5-11_ _S_ _4 _R1.fq.gz: 07 = individual, 5-11 = date, S = Sick, 1= stage1 (sampled on R1= right read bc they are paired reads) + R2 (left read)

zcat - program to open zipp files

``` 
zcat FILENAME | head
```


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

**Lab notebook tutorial:**

Typora: Command / = lets you togle back and forth from markdown to how it will look. 

To change Andrew's template use Typora, then use GitHub desktop and commit changes to master, then sync. 

____ send links to Andrew by Wednesday____ 



~/ = to go to my directory

cd scripts

Shows what is in the folder: 

```R
total 1932280

-rw-r--r--. 1 lcaicedo users 1926280714 Feb 13 11:14 075-11S4bwaaln.sam

-rw-r--r--. 1 lcaicedo users   25887148 Feb 13 11:11 075-11S4R1.fq.gz_left_clean_paired.fq.sai

-rw-r--r--. 1 lcaicedo users   26467948 Feb 13 11:13 075-11S4R2.fq.gz_right_clean_paired.fq.sai

-rwxr-xr-x. 1 lcaicedo users        903 Feb 13 11:07 bwaaln_lc.sh

-rw-r--r--. 1 lcaicedo users        835 Feb  8 10:44 cd

```



I have a SAM file 

To save tail: -n 100 07_05-11 _ S _ 4_bwaaln.sam > tail.sam

To open the tail: vim tail.sam

:set nowrap

```R
[lcaicedo@pbio381 scripts]$ tail -n 100 07_5-11_S_4_bwaaln.sam > tail.sam
[lcaicedo@pbio381 scripts]$ vim tail.sam

```

$ = to get to end in vim (where Tags are)

Tags to look for:

NM= edit distance

MD= mismatching positions/bases

X0= number of best hits

SAM files:



If seq don't map uniquely, they are discarded

```R
$ grep -c XT:A:U 07_5-11_S_4_bwaaln.sam 
648279
$ grep -c X0:i:1 07_5-11_S_4_bwaaln.sam 
650436
```

How many genes map to it? Use of python script

takes SAM file and counts how many reads uniquely mapped to reference sequence

Renaming of files

```R
sed -i 's/::/\_/g' 07_5-11_S_4_bwaaln.sam 
```

sed = 

's = search

/::/ = find

/\ _ / = replace

Primarily interested in UniqueTotalReads 

```R
#Go to the main directory and copy the script to my directory's scripts folder
cd /data/scripts
cp countxpression_pe.py ~/scripts      #or copy to your directory with the .sam file

#On my directory, inside scripts (w SAM file), run the python script
python countxpression_pe.py 20 35 countstatssummary.txt YOURFILENAME.sam
```

only pay attention to quality scores between 20 and 35

_____

02/22/2017

Navigate to data folder in terminal:

```R 
cd /data/project_data/DGE
ll
#shows 5 files in that file
#open up another terminal window, (navigate to the folder where I want to move it)
cd Documents/EcoGenomics/ 
scp lcaicedo@pbio381.uvm.edu:/data/project_data/DGE/* .
```

On Rstudio:

Install DEseq package:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
a
```

Then open script downloaded: DESeq2_exploreSSW_trim.R

___

02/27/2017: New Script for edited data

___

03/01/2017: Discussion of the models seen in last stcript

Line 225 of script: Is this interaction significant genome-wide

log(Likelihood (full model) - Likelihood (reduced))

​                           H,D,(not HxD)                  H,D

​	= LR ~ X^2, df = # of param full - # of param reduced

— Likelihood of the data given the model (sort of like Bayesian but without prior parameters)





-----

03/05/17

we have 94 reds, but 24 samples; we need to use one of two approaches towards analyzing each sample separately:

01_rep1.sam, 01_rep2.sam.. —> Merge using "samtools"

or

01_rep1.sam, 01_rep2.sam.. —> Call SNPs —> comare reps within individuals



vim head...

on vim- : set no wrap

applied 2 filters: any genotypes probablity called less than 95%

```R
vcftools --vcf SSW_bamlist.txt.vcf
```



*How could we quickly find out how many SNPs were flagged as unresolved?*

Filter for depth and quality

```R
grep "unres" SSW_bamlist.txt.vcf | wc
5631864 185851488 1028494934
grep "para" SSW_bamlist.txt.vcf | wc
   4354  143652  795592
```

wc= word count

Initial: 7.47 M

unresolved: 5.63 M (5631864)

parlogs: 4354

remaining: 1.8 M

— Anything not balletic is an error

```R
vcftools --vcf SSW_bamlist.txt.vcf --min-alleles 2 --max-alleles 2
--vcf SSW_bamlist.txt.vcf --maf 0.02
vcftools --vcf SSW_bamlist.txt.vcf --max-missing 0.8
```

allow 20% missing data

So far we have selected filters independently, 

place together:

enter r through terminal: R

```R
> getwd()
[1] "/data/users/l/c/lcaicedo"
> hardy <- read.table("out.hwe", header=T)
> str(hardy)
'data.frame':	442 obs. of  8 variables:
 $ CHR               : Factor w/ 111 levels "TRINITY_DN35598_c0_g1_TRINITY_DN35598_c0_g1_i1_g.5802_m.5802",..: 65 65 100 100 100 100 100 100 88 88 ...
 $ POS               : int  4566 4665 978 1404 1722 3426 3729 3912 115 141 ...
 $ OBS.HOM1.HET.HOM2.: Factor w/ 27 levels "10/11/3","11/0/13",..: 27 22 27 27 20 27 22 18 18 27 ...
 $ E.HOM1.HET.HOM2.  : Factor w/ 16 levels "10.01/10.98/3.01",..: 14 12 14 14 11 14 12 10 10 14 ...
 $ ChiSq_HWE         : num  0.0109 0.1067 0.0109 0.0109 0.1983 ...
 $ P_HWE             : num  1 1 1 1 1 1 1 1 1 1 ...
 $ P_HET_DEFICIT     : num  1 1 1 1 1 1 1 1 1 1 ...
 $ P_HET_EXCESS      : num  1 0.936 1 1 0.874 ...
> hardy[which(hardy$P_HET_EXCESS<0.001),]
> quit()
```

ID just loci with significant variation. 

______

03/08/2017

- login to server

  ```terminal
  cd /data/project_data/snps/reads2snps/
   32235925 Mar  7 22:21 SSW_byind.txt.vcf.gz
  ```

  The file is a gz file, compressed file (shown in red)

  ```
  vcftools --gzvcf SSW_byind.txt.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.02 --max-missing 0.8 --recode --out ~/SSW_all_biallelic.MAF0.02.Miss0.8
  ```

  —gzvcf (i will give you a zipped version of the vcf file)

  Edits: biallelic filter, biallelic freq. and eliminate missing data

  output saved to home directory (~) and given the name: SSW_all_biallelic.MAF0.02.Miss0.8 

  ```
  gzip SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf
  vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --hardy
  ```

  output in uncompressed format, changed to compressed

  vcf tools used for HW

  ``` R
  hwe <- read.table("out.hwe", header=TRUE)
  > str(hwe)
  #'data.frame':	1464 obs. of  8 variables:
  # $ CHR               : Factor w/ 363 levels #"TRINITY_DN27892_c0_g1_TRINITY_DN27892_c0_g1_i1_g.3123_m.3123#",..: 358 358 358 358 358 358 358 358 358 358 ...
  # $ POS               : int  4733 5850 5865 5869 5874 6096 #6146 6201 6289 6325 ...
  # $ OBS.HOM1.HET.HOM2.: Factor w/ 52 levels #"10/0/12","10/10/2",..: 36 51 36 36 36 36 36 29 36 48 ...
  # $ E.HOM1.HET.HOM2.  : Factor w/ 22 levels #"10.23/9.55/2.23",..: 14 1 14 14 14 14 14 11 14 1 ...
  # $ ChiSq_HWE         : num  0.0119 1.4547 0.0119 0.0119 #0.0119 ...
  # $ P_HWE             : num  1 0.361 1 1 1 ...
  # $ P_HET_DEFICIT     : num  1 0.954 1 1 1 ...
  # $ P_HET_EXCESS      : num  1 0.276 1 1 1 ...
  > summary(hwe)
                                                               CHR      
   TRINITY_DN45147_c0_g1_TRINITY_DN45147_c0_g1_i3_g.18680_m.18680:  34  
   TRINITY_DN46382_c0_g1_TRINITY_DN46382_c0_g1_i1_g.22149_m.22149:  28  
   TRINITY_DN45750_c0_g1_TRINITY_DN45750_c0_g1_i2_g.20209_m.20209:  27  
   TRINITY_DN47302_c3_g1_TRINITY_DN47302_c3_g1_i2_g.25471_m.25471:  21  
   TRINITY_DN46789_c1_g3_TRINITY_DN46789_c1_g3_i1_g.23393_m.23393:  20  
   TRINITY_DN46938_c1_g1_TRINITY_DN46938_c1_g1_i1_g.24007_m.24007:  19  
   (Other)                                                       :1315  
        POS         OBS.HOM1.HET.HOM2.        E.HOM1.HET.HOM2.
   Min.   :   1.0   21/1/0 :822        21.01/0.98/0.01:822    
   1st Qu.: 179.0   20/2/0 :202        20.05/1.91/0.05:240    
   Median : 321.0   19/3/0 : 96        19.10/2.80/0.10:103    
   Mean   : 630.5   18/4/0 : 69        18.18/3.64/0.18: 82    
   3rd Qu.: 728.2   17/5/0 : 44        17.28/4.43/0.28: 60    
   Max.   :6511.0   21/0/1 : 38        16.41/5.18/0.41: 29    
                    (Other):193        (Other)        :128    
     ChiSq_HWE             P_HWE           P_HET_DEFICIT      
   Min.   : 0.000094   Min.   :0.0000004   Min.   :0.0000004  
   1st Qu.: 0.011898   1st Qu.:1.0000000   1st Qu.:1.0000000  
   Median : 0.011898   Median :1.0000000   Median :1.0000000  
   Mean   : 0.943981   Mean   :0.9194100   Mean   :0.9216362  
   3rd Qu.: 0.117787   3rd Qu.:1.0000000   3rd Qu.:1.0000000  
   Max.   :22.000000   Max.   :1.0000000   Max.   :1.0000000  
                                                              
    P_HET_EXCESS      
   Min.   :0.0005731  
   1st Qu.:0.9767442  
   Median :1.0000000  
   Mean   :0.9432001  
   3rd Qu.:1.0000000  
   Max.   :1.0000000

  which(hwe$P_HET_DEFICIT<0.01) #rows with 0.01
   hwe[which(hwe$P_HET_DEFICIT<0.01),] #just these rows but all of the columns
  #Output SNPS in these transcripts(CHR)
  quit()
  ```

Calculate allele frequencies on vcftools

- by groups: Sick (SS, HS) and Healthy (HH)

```unix
cd /data/project_data/snps/reads2snps
vim ssw_healthloc.txt
#create 2 output files: for sick and healthy
grep "SS" ssw_healthloc.txt > ~/S_OneSampPerInd.txt
grep "HS" ssw_healthloc.txt > ~/S_OneSampPerInd.txt
grep "HH" ssw_healthloc.txt > ~/H_OneSampPerInd.txt
cat H_OneSampPerInd.txt 
[lcaicedo@pbio381 ~]$ cat H_OneSampPerInd.txt 
10	HH	INT	N
24	HH	INT	Y
27	HH	INT	Y
31	HH	SUB	Y
32	HH	SUB	Y
33	HH	SUB	Y
34	HH	SUB	N
35	HH	SUB	Y
[lcaicedo@pbio381 ~]$ cut -f 1 H_OneSampPerInd.txt > H_OneSampPerInd2.txt 
[lcaicedo@pbio381 ~]$ cat H_OneSampPerInd2.txt 
10
24
27
31
32
33
34
35
# use vcftools to calc allele freqs for H_OneSample...
[lcaicedo@pbio381 ~]$ vcftools --gzvcf SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz --freq2 --keep H_OneSampPerInd2.txt --out H_Allelefreq
# edit headers
vim H_Allelefreq.frq 
i #to edit
{FREQ} changed to MAJOR	MINOR
(esc) key, :wq
```

Use R

