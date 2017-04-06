

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

we have 94 reads, but 24 samples; we need to use one of two approaches towards analyzing each sample separately:

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

  —gzvcf (will give you a zipped version of the vcf file)

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

----

03/20/17

n=24

1. Final VCF data —> filter —> Output H.D.
2. Estimate allele frequency of H & S —> f(H) - f(S)
3. Fst between H vs S —> output to laptop —> plot on R
4. Estimate Pi(syn), Pi(nonsyn), ratio syn/nonsyn —> output to laptop, compare to Rominguier

wc = word count, counts number of rows and number of characters

```
$ grep "HS\|SS" ssw_healthloc.txt | cut -f1 >~/S_SampleIDs.txt
# "HS\|SS" match either HS or SS
```

---

03/21/2017

terminal: 

```
cd /data/project_data/snps/reads2snps
screen
/data/popgen/dNdSpiNpiS_1.0 -alignment_file=SSW_by24inds.txt.fas -ingroup=sp -out=~/dNdSpiNpiS_output
# keys: control+a+d
```

```
[lcaicedo@pbio381 reads2snps]$ cat SSW_bamlist.txt.sum
################################################################################
#                              Biological Summary                              #
################################################################################

Selected ingroup species: sp

Number of analyzed individual: 24 (from 1 population(s))

Total number of contig used for sequence analysis: 1113

Total number of SNPs: 5040

  - Biallelic: 4991
  - Triallelic: 49
  - Quadriallelic: 0

Fit: # like Fst

Average Fit: -0.0507419 [-0.06817; -0.031933]
(Fit calculated in 902 contigs)
# Hi Fit: inbreeding, high homozygocity

Weir & Cockerham Fit (Evolution 1984):

Average Weir & Cockerham Fit: 0.00703754 [-0.017669; 0.032047]
(Fit calculated in 902 contigs)

piN/piS ratio:

Average piS in focal species: 0.00585312 [0.005172; 0.006598]
# purifying selection, 
Average piN in focal species: 0.00154546 [0.00133; 0.001782]
Average piN / average piS: 0.264041 [0.223914; 0.310575]
#higher value means that selection isn't as effective at eliminating deleterious mutations.tends to be low for bacteria. nd invertebrates with high pop. sizes. to know if high or low we need to compare to Rominguier data(metazoa) a
(piS and piN calculated in 902 contigs of average length 50)

Individual heterozygosity:

H(03_5-08_S_2): 0.00292949
H(07_5-08_S_1): 0.00212127
H(08_5-08_H_0): 0.00160917
H(09_5-08_H_0): 0.00146331
H(10_5-08_H_0): 0.00118121
H(14_5-08_S_2): 0.00240668
H(15_5-08_H_0): 0.00153744
H(19_5-11_H_0): 0.00209038
H(20_5-08_H_0): 0.000972337
H(22_5-08_S_1): 0.00127272
H(23_5-17_S_2): 0.00252627
H(24_5-08_H_0): 0.000786939
H(26_5-08_S_2): 0.00145689
H(27_5-08_H_0): 0.00185711
H(28_5-08_S_1): 0.00137447
H(29_5-08_S_2): 0.00200665
H(31_6-12_H_0): 0.00178091
H(32_6-12_H_0): 0.00216789
H(33_6-12_H_0): 0.00210024
H(34_6-12_H_0): 0.00186592
H(35_6-12_H_0): 0.00229033
H(36_6-12_S_1): 0.00264985
H(37_6-12_H_0): 0.00230314
H(38_6-12_H_0): 0.0028128
```

Pi (s) -> estimator of theta 0 = 4Neu 

Ne = Pi(s)/ 4u

if: u~ 4X10^-9 (literature: Mike Lynch mutation rates of invertebrates)

Ne= 0.00555/ (4x 4X10^-9)= 365,625

-assumes generation time for sea stars is the same than estimated (lit)

- if gen. times are higher, then Ne wil be higher


---

Terminal: 

Have files on my directory: 

- SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.gz 
- ssw_healthloc.txt

Fetch to my computer, use script on tutorial: https://adnguyen.github.io/2017_Ecological_Genomics/Tutorial/2017-03-22_Tutorials_PopulGenomics4.html

-----

03/29/17: lab tutorial from ADMIXTURE (https://adnguyen.github.io/2017_Ecological_Genomics/Tutorial/2017-03-22_Tutorials_PopulGenomics4.html)

dataset of i individuals @ j SNPs -> Pr(G|K,Q,P)

- Q= ancestry coefficient (proportion of ancestry in each population in an individual genome) 
- P= allele frequencies at each of the K populations

ADMIXTURE: analysis is not Bayesian

- to choose a model:
  - If G= 24 individuals: [ matrix divided in chunks and masks one of these] Cross validation: If I develop a model on all except the masked one, how good is the model at predicting the genotype of the masked one.(similar to a Jacknife), we want a low cross validation score, choose the K that generated the best cross validation. 
- If A= reference allele, T = alternate allele

Import

copy to home directory the files: vcf2admixture_SSW.spid, SSW_tidal.pops, vcf2geno.sh, and ADMIX.sh

```
vim vcf2geno.sh
# change the input file name (SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf)and output file name (SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno)
./vcf2geno.sh # to bash the file
vi ADMIX.sh #vi same as vim
--
#!/bin/bash

# Run ADMIXTURE to determine the number of genetic clusters in the SNP data, 
# and the ancestry proportions of each individual

# Here's the utility of 'for loops'...

for K in {1..10}

do

admixture -C 0.000001 --cv #C: is small, the smaller the more you can make the algorithm climb to find the peak# ./SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno $K \
| tee log${K}.out # |= take the output from admixture analysis and send to a log file into .out file, tee= take everything#

done
# $ = 
# After the for loop finishes, you can use 'grep' to grab the values of the CV from each separate log file and append them into a new summary text file.

grep CV log*.out >chooseK.txt # CV= cross validation error and send this bit to new file called chooseL.txt#
```



write and exit vim

bash ADMIX.sh

- to see output Ks: cat chooseK.txt:

```
log10.out:CV error (K=10): 0.54329
log1.out:CV error (K=1): 0.52967 #this is the lowest value, so one population#
log2.out:CV error (K=2): 0.62210
log3.out:CV error (K=3): 0.72173
log4.out:CV error (K=4): 0.81377
log5.out:CV error (K=5): 0.89176
log6.out:CV error (K=6): 0.94948
log7.out:CV error (K=7): 0.89471
log8.out:CV error (K=8): 0.96684
log9.out:CV error (K=9): 0.95880
```

If, K:1-20: gives lowest value for K=20 because of overfitting

--------

HW: take tools and try with dif. filtering strategies: 

```
--min-alleles 2
--max-alleles 2
--max-missing 0.8, or 1.0 (eliminate all missing data)
--recode
--out filename
--maf: minor allele frequencies (seen in Gompert) can go as low as 0.01, ussually 0.05
--minDP: read depth (low being 5)
--hwe: based on P-value, eg. 0.01 elliminate any that violate HWE, to explore effect of HW

analyze 2 vcf files using PCA, DAPC or ADMIXTURE (choose one)

# to remove individuals:
vcftools --vcf filename.vcf --remove-indivs 03 --remove-indivs 17 --min_alleles 2 --max-alleles 2
# OR have a text file with the individuals to remove:
--remove list_of_indiv_to_remove.txt ...
# how to select individuals to remove? indiv. above a threshold of missing data (eg. 10%) OR using a heat map with

```

* use journals from class or look at Molecular ecology.

_______

04/03/17

search in vim: /whatIwanttosearch



----

**Why might we want to BLAST to different databases?**

more prob of hit from NR, but get more functional information with Uniprot

​	<u>Query</u>			    <u>Subject</u>

​	.cds -> DIAMOND      nr

ORF .pep -> blastp	    uniprot -> uniprot ID-> GO IDs KEGG taxa

​	.pep-> blastp		    P.miniata.pep & S.purpuratus.pep->