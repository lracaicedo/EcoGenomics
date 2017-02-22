
# 2017 Ecological Genomics Course

### Author: Laura Caicedo-Quiroga     


## Overall Description of notebook      
Here you can see my online notebook for the Ecological Genomics course offered for the Spring of 2017. I will take notes on themes on the subject, paper discussions and methods discussed on the lab portion of the class. 

## Date started: (2017-01-18)   
## Date end:   -    




<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.  


### Table of contents        
* [Page 1: 2017-01-18](#id-section1). Introduction to Ecological Genomics
* [Page 2: 2017-01-20](#id-section2). Day 2 (Info update rubric, course notes, paper discussion)
* [Page 3: 2017-01-25](#id-section3). Day 3
* [Page 4: 2017-01-30](#id-section4) . Day 4: Group projects
* [Page 5: 2017-02-01](#id-section5). Week 3, Day 4, class notes , Group presentations of project ideas
* [Page 6: 2017-02-01](#id-section6). Week 3, Day 5, command line stuff
* [Page 7: 2017-02-03](#id-section7). Installing trinity onto UVM cluster
* [Page 8: 2017-02-06](#id-section8). Week 4, Day 6, RNA-seq
* [Page 9: 2017-02-08](#id-section9). Week 4, Day 7, RNA-seq cont'd + paper discussion DePanis et al. 2016; MolEco
* [Page 10: 2017-02-10](#id-section10). Prepping for leading journal club discussion: 2015-02-15; Zhao et al. 2016; *MBE*
* [Page 11: 2017-02-13](#id-section11). Week 5, Day 8, RNA-seq mapping and paper discussion: Johnston et al. 2016, *Molecular Ecology*
* [Page 12:](#id-section12).
* [Page 13:](#id-section13).
* [Page 14:](#id-section14).
* [Page 15:](#id-section15).
* [Page 16:](#id-section16).
* [Page 17:](#id-section17).
* [Page 18:](#id-section18).
* [Page 19:](#id-section19).
* [Page 20:](#id-section20).
* [Page 21:](#id-section21).
* [Page 22:](#id-section22).
* [Page 23:](#id-section23).
* [Page 24:](#id-section24).
* [Page 25:](#id-section25).
* [Page 26:](#id-section26).
* [Page 27:](#id-section27).
* [Page 28:](#id-section28).
* [Page 29:](#id-section29).
* [Page 30:](#id-section30).
* [Page 31:](#id-section31).
* [Page 32:](#id-section32).
* [Page 33:](#id-section33).
* [Page 34:](#id-section34).
* [Page 35:](#id-section35).
* [Page 36:](#id-section36).
* [Page 37:](#id-section37).
* [Page 38:](#id-section38).
* [Page 39:](#id-section39).
* [Page 40:](#id-section40).
* [Page 41:](#id-section41).
* [Page 42:](#id-section42).
* [Page 43:](#id-section43).
* [Page 44:](#id-section44).
* [Page 45:](#id-section45).
* [Page 46:](#id-section46).
* [Page 47:](#id-section47).
* [Page 48:](#id-section48).
* [Page 49:](#id-section49).
* [Page 50:](#id-section50).
* [Page 51:](#id-section51).
* [Page 52:](#id-section52).
* [Page 53:](#id-section53).
* [Page 54:](#id-section54).
* [Page 55:](#id-section55).
* [Page 56:](#id-section56).
* [Page 57:](#id-section57).
* [Page 58:](#id-section58).
* [Page 59:](#id-section59).
* [Page 60:](#id-section60).

------
<div id='id-section1'/>
### Page 1: 2016-01-18. Ecological genomics: Introduction

#### **Uses**    
- Understand the genetic source of adaptations
- Approaches and tools for studying relationship: <u>genome/phenotype/environment</u>
  - Ecological issues such as nutrient cycling, population structure, life history vairation , trophic interaction, stress responess, and adpatation to environmental change 
- Environmental change and adaptation studies

* On any system, not only model organisms 
* Next-gen sequencing: increasingly useful but with huge datasets (computational challenge)    

#### Questions

- Which genes are expressed at certain times?
- What is the influence of the microbiome?
- What can we find on environmental samples through metagenomes?
- What is the influence of epigenetics on evolution? (independence from genetic code)
- What are the patterns of genetic diversity for Popgen?
- What are the constraints of genomes?

#### Methods

- De-novo genome assembly
- Differential expression analysis: RNA-Seq
- Microbial community diversity, 16S metagenomic sequencing
- If reduced genome or one sample, RAD-Seq

#### Proccesses  studied

* Speciation, hybridization
* Local adaptations        

##### Melissa's and Steve's stories [see Andrew's notes] (https://github.com/adnguyen/Notebooks_and_Protocols/blob/master/2017_Eco_Gen_ANBE_nb.md#id-section1)  

------
<div id='id-section2'/>

### Page 2: 2017-01-23. Day 2, course notes

#### Info update rubric

- Outline : can give handouts
- 20 minutes
- learning/ engaging activity
- use of board effectively
- take home messages (star)
- examples from literature (beyond those given onBb)
- Glossary is built here

#### <u>Melissa's Info update</u>

##### Outline

1. Advances in Seq. Technologies
2. Range of Applications
3. General Library Prep. Workflow
4. Sequencing-by-synthesis (SBS)
5. Other Technologies
6. Learning Activity



1) Sequencing Technologies

- Human Genome Project: 2001-2003 ABI – Sanger, took 15 yearsand sequenced 1 genome, $3B cost
- Development of Hi Seq Tten in 2014 by **Ilumina**. This tool 1 day and 45 whole human genomes and cost $1000each. Process: Microarray size slide: Ilumina flow cell with 8 lanes werecapillaries are sampled by billions. 

2) Applications: **90%global data**, mostly Ilumina seq. 

- WGS (whole genome sequencing), RNAseq, ChIPseq (Chromatinimmunoprecipitation, protein bound transcription factors), targeted/capturesequence (use of probe= stretch of DNA that you synthesize, therefore you needto know about the target)

3) Workflow

- Kind of library prep to make:
  - Where is genetic variation? (Phenotypes)
  - number of samples: population, individual,comparative studies
  - Model or not
  - Demographic history of population
  - Adaptive genetic variation
  - Gene expression variation
  - Determined by:
    - length ofreads: longer reads will be easier to assemble 
    - number of reads
    - distribution of the reads: knowing relativeposition helps with assembly, are they randomly distributed or not
- Steps:
  - Extraction: DNA or RNA -> **cDNA**
  - Fragment sample
  - Ligate adaptors: individual barcode
  - Add sequence adaptors
  - PCR amplification 
- Reduced representation
  - RNA: coding regions
  - GBS/RADseq:near restriction sites

 4) SBS (Sequencing by Synthesis)

-      Flow cell ,has lanes and oligos (eg. P7, P5) that attaches to your DNA that has adaptors thatsample by ligation. Your DNA is attached to a sequence adaptor and may have anID barcode. Then **bridge amplification**where the pieces bend over and amplify (to copy and make signal stronger). Bridgereleased and then cluster generation (clusters of the same sequence). For eachcluster a snapshot is taken of the nucleotides, as the polymerase is adding thenucleotides. 

-      PacBio: Single Molecule in Real Time. Givesreally long reads but accuracy is lower. 

* Drawings on notebook*

6) Learning Activity

### Glossary:

Ilumina reads: short (50bp), long (100,150,300 bp), 

SMRT: extra long (10000 – 60000 bp). Can also be single vs.paired end. 

______

#### <u>Paper discussion:</u> Genome sequencing and population genomics in non-model organisms (Hans Ellegren, 2013)

- Omics: when is it useful? Can be used as an unbiased view
- Storage for genomic data! Since it should be publicly accessible (Short read archive) - NCBI
- Quality control: files are raw, quality processing is done afterwards
- Table 1: 2,200 update of reference genomes of eukariotes (600 on paper)
- Great proportion of non-model organisms. Does not include transcriptomes. 
- How do we determine what is a reference genome?  In scaffolds, not in chromosomes. There is a database to query. Sample size is 1 (not always) but how to know if it is a reference that works for a specific question?

 — Next Monday: 4 Info-updates— 

------
<div id='id-section3'/>

### Page 3: 2017-01-25. Day 3, course notes    

#### <u>Info update: QTNs</u>

**Outline**

1. What are QTN??

   1) QTN: "Quantitative Trait Nucleotides"

   ```Examples:
   eg.
   - flowering time: Continuos/ Quantitative traits
   - flower colors: discreet/ Mendelian traits
   - thermal tolerance
   - venom potency
   - altitude tolerance
   - defense compounds
   - toxin tolerance
   - draught tolerance
   - altitude tolerance

   ```

   - Quantitative and continuous
   - Quantitative traits: you can decompose in genotypes, each genotype contributes an ammount to the phenotype. 


- Fisher: Small effect vs. long effect mutations (see notebook)

2. Quantitative Genetic Theory of Adaptive Traits

   - Va
   - h^2

3. Methods

   - QTL - Linkage mapping (forward genetics)
   - GWAS- Genome-wide association studies: "fancy regression" (forward genetics)
   - Selection scans (reverse genetics)

#### Paper discussion:

<u>Rockman, (directed by Allison Brody)</u>

- are QTNs useful?

3 points presented: 

1. small effect are the vast majority (mendelian), while the allelic frequency we are drawing from are not the most representative. "unmeasurable interesting traits"
2. Theory
3. small effect alelles do not operate in the same manner: bias in effects. 

impacts of the large effect and small effects

<u>Lee et al. 2014 (directed by Lisa Chamberland)</u> 

- Still useful, we might get to the small-effect QTNs with these methods

Discussion:

Use in combination with other methods

Usefulness in id. If QTN's for pop structure, and the genetic basis of an adaptation that may not apply to different populations. 

Is it worth it for the scope of a PhD program?

Communication between findings and researchers that will use this knowledge.



#### Student project introduction: Melissa 

—notes on data paper handed— 

------
<div id='id-section4'/>

### Page 4: 2017-01-30. Week 3, Day 4

##### <u>Group projects</u>

- Intertidal vs Subtidal: Subtidal much easier to collect because they are fully submerged. Therefore gene expression and 16S (microbiome) differences because of handling effect (could be stress responses), not using snp data. One month time lag between collections, May- June.
- Sea-star associated Densovirus- single stranded dna virus. QPCR for all samples: negative for densovirus. Samples from epidermal tissue
- sp. winter spawning
- Sea-star biology: only identification juvenile vs adult not age class. Sex id. with dissection
- RNA-seq: anything with a Poly-A tail will be amplified. Could be fungi or annelids...

#### Group B: Intertidal vs. Subtidal

- gene expression and microbiome differences

  ​

Gene expression: Differences in expression between the groups (identify a treshold of (significant) expression to determine high or low expression) and then pinpoint genes of interest. Congruence and difference in the genes expressed. 

Level of expression of all genes by groups.Differential expression. What about house keeping genes?

Microbiome differences: ribosomal RNA , crude estimates of abundance and diversity of microbiota. <u>Community structure</u> of the micro biome (multid. scaling, or discriminant analysis).  Diversity: would decrease with illness. Hypotheses: more diverse microbiota leads to healty ind. or specific community composition towards healthy individuals. 

Samples to use: Day 3 only

Next steps: 

Lit search on Microbiome communities on sea star epidermis, wasting disease, inmune response on sea stars.

— extra notes on Handout — 

*Write up ideas into 1 page proposal using guideline provided and integrating feedback from the discussion. Include specific samples. For Monday, MS-Word and e-mail to Steve and Melissa*



Wednesday: 4 library prep types Blitz

------
<div id='id-section5'/>
### Page 4: 2017-02-17. Day 5

#### Announcements

- Send link to github online notebook to adnguyen@uvm.edu, 2 "notebooks" notes and coding
- Sign-ups
- Project proposal due next Monday by email: one per group
- Transcriptomics next week



#### <u>WGS</u>: Info update by me {notes on notebook}

- everything

#### <u>RNA-seq</u>: Info update 

#### gene space expressed

1. Advantage
   1. Differential gene expression: in tissues there is variation
   2. Allele specific expression: environmental response of adaptation
   3. Functional relevant subset of the genome
   4. vs. Micro-array:  RNA seq has wide range of expression value, not in micro-array. No saturation of analog-type signal (flourescent) in micro-array . Informs on spilicing events. 
2. Limitation
   1. Post transcriptional events not identified
3. Work-flow
   1. Set-up
   2. Wet-lab
   3. Sequencing strategy
   4. Bio-informatics
   5. Statistical Measures

Work Flow: 

1. Set-up

   1. Protein coding or regulatory non-coding
   2. Reference genome present?
   3. Alternative splicing?
   4. Technology for research question
   5. Population or specific treatment?
   6. Stat: Biological Replication
   7. Choice of Tissue: circadian rythim. If small organism: use pooled samples

2. Wet-Lab

   1. Sample: RNA extraction for good yield

      1. RNase free environment
      2. Treat with DNase (to avoid DNA hybrids)
      3. Get rid of rRNAs

   2. cDNA: 5'————AAA

      ​            3'————TTT oligo dT primer

      ​	Reverse Transcriptase

      ​	3'———TTT-5'		Library: single end and paired end

      ​	5'———AAA-3'  (cDNA)

   3. Platform: Pyrosequencing by Roche, Ion Torrent. GA?Hiseq by Ilumina (GC content)

      1. Error Profiles of the platforms: Incorrect homopolymers

      2. Sequence coverage: >100 millian bp

      3. Programming

         ​

      Experimental Set-up -> Tissue Prep + Library -> High throughput seq. -> Transcriptome reconstruction -> Alignment of Reads: - Read qualification - Biological Influence -Marker development

de novo sequencing: 

. . .  . 

5'(UTR) CDS___3' UTR

Assembly: -^-^-^ Reconstruction and slice junction



#### <u>Amplicon Sequencing</u>: one gene

Def: Targeted approach for analyzing genetic variation in specific genomic regions (Amplicon= targeted gene region to be amplified via PCR w/ specific primers)

- Methods

  - Libray Prep: 

    - Extract DNA (specific region)
    - PCR 1 (16S) rRNA: Amplify gene, specific primers. —> Clean (use a gel) and can extract through the gel or use column and sequence if ok if not: —> PCR 2: barcodes and adaptors —> Clean again —> Pool —> Sequence

  - Sequencing

    - Platforms: 454, MiSeq (Ilumina)

  - Data Analysis (for 16S)

    - Trim adaptor sequence (e.g.. AGGTTT 7bp)
    - Align Data: use conserved areas, tell your program the seq to use 

    ​

- Applications

  - Id species
  - one gene of interest
  - Benefits: if you know the gene you want, you save time and money,

- Limitations: 

  - Computational

#### <u>GBS</u>: 

- Genotyping by sequencing at the same time (before this you didn't have to develop primers and then sequence). between RNA-seq and Amplicons

"Reduced Representation"

<u>RAD seq</u>: Resticton enzyme assisted DNA sequencing

- Lots of individuals from 
- lots of SNPs across the genome
- Don't care about specific genes



------
<div id='id-section6'/>
### Page 5: 2017-02-06. Day 6 

#### <u>Info update: RNA-seq</u> 






------
<div id='id-section7'/>
### Page 7: 2017-02-03. Installing trinity into the cluster 

Get the tar and zip file here: https://github.com/trinityrnaseq/trinityrnaseq/releases; 

[Installation instructions](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing%20Trinity)



1. Working directory: /data/popgen
2. Used "wget"
   * Downloading the file: 

​```
			wget "https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.3.2.tar.gz"
​```

3. Gunzip and tar commands to get executable.

​```
gunzip Trinity-v2.3.2.tar.gz
tar -xvf Trinity-v2.3.2.tar
​```

4. cd'd into triniy directory and use "make"

5. make plugins

6. Check to see if it works

   ```
   cd sample_data/test_Trinity_Assembly/
   ./runMe.sh
   ```




------



<div id='id-section8'/>
### Page 8: 2017-02-06. Week 4, Day 6, RNA-seq

## Class outline:

1. info update with Melissa on RNA-seq
2. Paper discussion
3. RNA-seq pipeline through the server



## 1 .info update with Melissa on RNA-seq

Pre-overall workflow:

1. Approach
2. Experimental design
3. Library Prep
4. Sequencing facility
5. Receive data
6. Computer/server setup 
7. Processing

==Processing Workflow:==



​```mermaid
graph TD
A{Raw reads, fast-q} --> B[Clean reads: Get rid of adaptors]
B --> C[Clean reads: Check nucleotide quality]
C --> D[Clean reads: Length? ]
D --> E(Evaluate Quality)
E --> F{De novo Transcriptome Assembly, fasta}
F --> G[Evaluate Transcriptome assembly]
G --> Y(compare to closely related species)
G --> P(CEG)
G --> W(N50 and number of contigs)
G --> y(BLAST to annotate)
F --> K(Map raw reads, SAM files)
A --> K
K --> H(Extract read count info, gene expression)
H --> U(DGE analysis)
K --> V(Identify SNPs)
V --> J(Population genomics analyses)
J --> a(Genetic Differentiation)
J --> b(Population Structure)
J --> h(test for signatures of selection)

​```



## 2. Paper discussion: Dunning et al. 2014

Title: Divergent transcriptional responses to low temperature among populations of alpine and lowland species of New Zealand stick insects (Micrarchus).



Shared response to cold shock? Stick insects in new zealand. They all responded to cold shock differently through differential gene expression (~ 2000 unique genes). 



**Hypothesis**: We hypothesize that species with poor dispersal ability are likely to have strong phylogeographic structure as a result of genetic drift and, possibly, local adaptation.

* The resulting divergent genetic backgrounds are likely to contribute to variation in the intra- and interspecific transcriptional responses to environmental stress.

They're poor dispersers, so that you'd expect for them to be locally adapted, particularly to the environment. 



Approach: 

* They used SNPs and COI(mitochondria) to measure extent of population/species differentiation.
* De novo assembly in transcriptome
* 4 taxa; generated transcriptome data from experiment, 



What is phred quality? Quality score on the propbability of the error.  >30  chance it is "right".



What are unigenes? *Trinity components containing clusters of ‘contigs’ representing splice variants of the same locus.*



## RNA-seq tutorial

1. CD to path

   * ```UNIX
      cd project_data/fastq/
     ```

2. Files: 

   * ```UNIX
     ls
     07_5-08_S_1_R1.fq.gz  27_5-11_H_0_R2.fq.gz
     07_5-08_S_1_R2.fq.gz  27_5-14_H_0_R1.fq.gz
     07_5-11_S_4_R1.fq.gz  27_5-14_H_0_R2.fq.gz
     07_5-11_S_4_R2.fq.gz  27_5-20_H_0_R1.fq.gz
     08_5-11_S_1_R1.fq.gz  27_5-20_H_0_R2.fq.gz
     08_5-11_S_1_R2.fq.gz  28_5-08_S_1_R1.fq.gz
     08_5-14_S_1_R1.fq.gz  28_5-08_S_1_R2.fq.gz
     08_5-14_S_1_R2.fq.gz  28_5-11_S_1_R1.fq.gz
     09_5-08_H_0_R1.fq.gz  28_5-11_S_1_R2.fq.gz
     09_5-08_H_0_R2.fq.gz  28_5-17_S_2_R1.fq.gz
     09_5-14_S_2_R1.fq.gz  28_5-17_S_2_R2.fq.gz
     09_5-14_S_2_R2.fq.gz  29_5-08_S_2_R1.fq.gz
     10_5-08_H_0_R1.fq.gz  29_5-08_S_2_R2.fq.gz
     10_5-08_H_0_R2.fq.gz  29_5-11_S_2_R1.fq.gz
     10_5-11_H_0_R1.fq.gz  29_5-11_S_2_R2.fq.gz
     10_5-11_H_0_R2.fq.gz  29_5-14_S_2_R1.fq.gz
     10_5-20_S_2_R1.fq.gz  29_5-14_S_2_R2.fq.gz
     10_5-20_S_2_R2.fq.gz  31_6-21_H_0_R1.fq.gz
     15_5-17_S_3_R1.fq.gz  31_6-21_H_0_R2.fq.gz
     15_5-17_S_3_R2.fq.gz  32_6-15_H_0_R1.fq.gz
     19_5-11_H_0_R1.fq.gz  32_6-15_H_0_R2.fq.gz
     19_5-11_H_0_R2.fq.gz  32_6-18_H_0_R1.fq.gz
     19_5-17_H_0_R1.fq.gz  32_6-18_H_0_R2.fq.gz
     19_5-17_H_0_R2.fq.gz  32_6-21_H_0_R1.fq.gz
     19_5-20_S_5_R1.fq.gz  32_6-21_H_0_R2.fq.gz
     19_5-20_S_5_R2.fq.gz  33_6-12_H_0_R1.fq.gz
     20_5-14_H_0_R1.fq.gz  33_6-12_H_0_R2.fq.gz
     20_5-14_H_0_R2.fq.gz  34_6-12_H_0_R1.fq.gz
     22_5-08_S_1_R1.fq.gz  34_6-12_H_0_R2.fq.gz
     22_5-08_S_1_R2.fq.gz  34_6-18_H_0_R1.fq.gz
     22_5-11_S_1_R1.fq.gz  34_6-18_H_0_R2.fq.gz
     22_5-11_S_1_R2.fq.gz  35_6-15_H_0_R1.fq.gz
     23_5-17_S_2_R1.fq.gz  35_6-15_H_0_R2.fq.gz
     23_5-17_S_2_R2.fq.gz  35_6-18_H_0_R1.fq.gz
     24_5-08_H_0_R1.fq.gz  35_6-18_H_0_R2.fq.gz
     24_5-08_H_0_R2.fq.gz  36_6-12_S_1_R1.fq.gz
     24_5-14_H_0_R1.fq.gz  36_6-12_S_1_R2.fq.gz
     24_5-14_H_0_R2.fq.gz  36_6-15_S_2_R1.fq.gz
     24_5-17_H_0_R1.fq.gz  36_6-15_S_2_R2.fq.gz
     24_5-17_H_0_R2.fq.gz  36_6-18_S_3_R1.fq.gz
     24_5-20_H_0_R1.fq.gz  36_6-18_S_3_R2.fq.gz
     24_5-20_H_0_R2.fq.gz  38_6-18_S_2_R1.fq.gz
     26_5-08_S_2_R1.fq.gz  38_6-18_S_2_R2.fq.gz
     26_5-08_S_2_R2.fq.gz  38_6-21_H_0_R1.fq.gz
     26_5-11_S_3_R1.fq.gz  38_6-21_H_0_R2.fq.gz
     26_5-11_S_3_R2.fq.gz  38_6-24_S_5_R1.fq.gz
     27_5-08_H_0_R1.fq.gz  38_6-24_S_5_R2_fastqc.html
     27_5-08_H_0_R2.fq.gz  38_6-24_S_5_R2_fastqc.zip
     27_5-11_H_0_R1.fq.gz  38_6-24_S_5_R2.fq.gz
     ```

   * Each class member will do R1 and R2(left and right reads)

   * number relates to scale of symptoms

   * Samples I'm doing!

   * ```UNIX 
     20_5-14_H_0_R1.fq.gz & 205-14_H_0_R2.fq.gz
     ```

3. look at our files: 

   * ```UNIX
     zcat 20_5-14_H_0_R1.fq.gz | head

     ```

   * Output:

     @J00160:63:HHHT2BBXX:1:1101:27823:1244 1:N:0:TCCGGAGA+AGGCTATA

     GNGCGTTATTATATGGTTTTATCTTCATTTNTTAAATGAACTTGATCTTGAATTTTTTTTTTTTTTTTTTTGGGGGATCGGAAGAGCACACGTNTGAACTC

     +

     A#AFFAFJJJJJJJJJJJJJJJJJJJJJJJ#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAAFFF-A<JFAJF<JJ-FJJJF#JF-J77A

     @J00160:63:HHHT2BBXX:1:1101:28635:1244 1:N:0:TCCGGAGA+AGGCTATA

     ANTGAGTAGAAGGAATCGGTCCACCATAAANAAGTGGAGGTTCCACATGGGCAAAGATGCCGGTACCATTCTTAACACTAGAAGAAGGAGCTTTTTCACTA

     +

     A#AFFJJJJJJJJJJJJJJJJJJJJJJJJJ#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJJJJJJJJJJJJJJJJJJJJJJJJAFJJJJJF

     @J00160:63:HHHT2BBXX:1:1101:29244:1244 1:N:0:TCCGGAGA+AGGCTATA

     GAGGCACTCATACAGGTTACACAGCTGAGANTAATTTATATCATATACTATAATGCATAATACATGTAAGCATCTCTATTGCTACATTGCCTGGTTATACA

4. This is a fastq file

   * Unique indentifier
   * ATGC's are sequence data for the reads; same length
   * 3rd line is a plus, just an indicator; before the quality score (PHRED) for each nucleotide bases (ASCII format to have double digit value to be represented by a single character). Makes downstream processing easier
   * Letters are good quality, non-letters are not.
   * greater than 30 is usually the cutoff

5. move into my directory and save the scripts into the scripts folder

   * ```UNIX
     cd scripts/
     ```

   * ```UNIX
     cp /data/scripts/trim_example.sh .

     ```

   * ​



## Actual script to run trimmomatic

​```
#!/bin/bash
      java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
                -threads 1 \
                -phred33 \
                 /data/project_data/fastq/20_5-14_H_0_R1.fq.gz \
                 /data/project_data/fastq/20_5-14_H_0_R2.fq.gz \
                 /data/project_data/fastq/cleanreads/"20_5-14_H_0_R1_clean_paired.fa" \
                 /data/project_data/fastq/cleanreads/"20_5-14_H_0_R1_clean_unpaired.fa" \
                 /data/project_data/fastq/cleanreads/"20_5-14_H_0_R2_clean_paired.fa" \
                 /data/project_data/fastq/cleanreads/"20_5-14_H_0_R2_clean_unpaired.fa" \
                 ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
                 LEADING:28 \
             TRAILING:28 \
             SLIDINGWINDOW:6:28 \
             HEADCROP:9 \
             MINLEN:35 \
​```





## Running script:



​```
./trim_example.sh &
​```



Output



​```
Input Read Pairs: 13876156 Both Surviving: 11211956 (80.80%) Forward Only Surviving: 2046413 (14.75%) Reverse Only Surviving: 257037 (1.85%) Dropped: 360750 (2.60%)
TrimmomaticPE: Completed successfully
​```





### Course notes and thoughts:

* Students may benefit from instructor led demonstration of implementing an electronic notebook.
  * It is also more engaging when the instructor takes notes or copy and paste code into a notebook
* Give more clear instructions for how journal club leaders should lead a paper. 
  * Students should give questions, rational, hypotheses, experimental design,
* ​

------
<div id='id-section9'/>
### Page 9: 2017-02-08. Week 4, Day 7, RNA-seq cont'd + paper discussion DePanis et al. 2016; MolEco

Transcriptomics take 2: 

## Info update by Lisa Chamberland

1. Intro
   * studying relationship between organism and environment
   * adding genomics you can look at this at a finer scale
   * observe rapid responses to the environment (gene expression)
   * advantages: can study wild systems (non-models; or non-traditional)
   * find silent genes that get expressed under novel stimuli
   * variation between individuals, within and among populations, 
   * novel transcripts without homologues in closely related organisms
   * Main methods:
     * microarrays
     * nexgen sequencing
2. Breif overview
3. Main questions
   1. How much variation is there in gene expression and how is it structured? (288 studies)
      * Heritable variation
      * Epigenetics
      * Wild populations
      * Qst-Fst comparisons
      * eQTL
      * Macroevolution
        * bottlenecks, drift, selection
   2. How do environmental stimuli affect gene expression? ( 136)
      * abiotic stress
      * environmental heterogeneity
      * host-parasite interactions
      * selective abiotic and biotic interactions
      * molecular basis of response
        * phenotypic plasticity
        * among genotypes
      * drawbacks : need to flash freeze and the transcriptome is just a snapshot
      * time course analyses
   3. how does gene expression affect phenotype? (15 studies)
      * alternative phenotypes
      * move from correlation to causation
        * transgenics, RNAi, crispr-cas9
4. future directions
   * combined microarrays and rna-seq 
   * Database for proposed ecologicla variations
   * Problems:
     * bias in signals 
     * ~~heterologous arrays~~
     * polyploidy
     * RNA pooling
     * statistical anlayses



Glossary: 



## Paper discussion

Background: Desert drosophila that grows in arid conditions. It has 2 hosts. The two cactus differ in alkyloids. 

Q: How does environment influence gene expression? 





## Coding session

Learning goals (skills):

1. scripts
   * executing scripts
2. paths
   * program
   * input
   * output
   * filenames in and out
3. moving through the directories
4. moving files
   * on the server 
   * ==from server to our computer (scp)==
5. ​

Practical:

1. Finish up our cleaning: trimmomatic
2. fastqc (visualizing)
3. make a table of # reads
4. design assembly tests
5. Start assemblies
6. Evaluate assembly



<u>Moved the trimmomatic files into this path:</u>



​```
/users/a/d/adnguyen/mydata/2017-02-08_cleanread
​```

<u>Ok now, run fastqc to check files</u>



​```
fastqc 20*
​```



<u>What do I get?</u>

​```
[adnguyen@pbio381 2017-02-08_cleanreads]$ ls
20_5-14_H_0_R1_clean_paired.fa
20_5-14_H_0_R1_clean_paired.fa_fastqc.html
20_5-14_H_0_R1_clean_paired.fa_fastqc.zip
20_5-14_H_0_R1_clean_unpaired.fa
20_5-14_H_0_R1_clean_unpaired.fa_fastqc.html
20_5-14_H_0_R1_clean_unpaired.fa_fastqc.zip
20_5-14_H_0_R2_clean_paired.fa
20_5-14_H_0_R2_clean_paired.fa_fastqc.html
20_5-14_H_0_R2_clean_paired.fa_fastqc.zip
20_5-14_H_0_R2_clean_unpaired.fa
20_5-14_H_0_R2_clean_unpaired.fa_fastqc.html
20_5-14_H_0_R2_clean_unpaired.fa_fastqc.zip
​```



<u>Now, I'll move it back to my computer: So in my home computer, this is my working directory.</u>



​```
/Users/andrewnguyen
​```



<u>Now we can move files from server to my home computer</u>

​```
scp adnguyen@pbio381.uvm.edu:~/mydata/2017-02-08_cleanreads/*.html .

​```

<u>This is what i get</u>

​```
andrewnguyen$ ls
20_5-14_H_0_R1_clean_paired.fa_fastqc.html
20_5-14_H_0_R1_clean_unpaired.fa_fastqc.html
20_5-14_H_0_R2_clean_paired.fa_fastqc.html
20_5-14_H_0_R2_clean_unpaired.fa_fastqc.html
Applications
Desktop
Documents
Downloads
Dropbox
Google Drive
Library
Movies
Music
Pictures
Public
Sites
Teaching
bower_components
node_modules
package.json
zScience

​```



**When doing an assembly. You can concatenate reads and then run fastqc.**



## Now we can run an assembly:

What can influence it?  

* number of individuals? 





instructor meeting:

* Coming up is assembly and mapping sequence reads. 


* There is a catch up day. We're ~ half a day behind. 
* Melissa make better quality assembly
  * Evaluate the assembly, 
  * Start annotation mapping (can use Trinitate; built into Trinity)
* Melissa (for Monday): finish cleaning all the samples and then make a good reference
  * Health and sick individuals; 50 million reads
  * BWA-mem; short mapping ==(by monday)==
    * students will map their own files —> sam files
    * then concatenate sam files
* for transcriptomics, they can extract read counts ffrom their sam files ==(Wednesday)==. (this involves a python script)
  * DEseq2 in R for transcriptomic analyses


------
<div id='id-section10'/>
### Page 10: 2017-02-10. Prepping for leading journal club discussion: 2015-02-15; Zhao et al. 2016; *MBE*

reference:

Zhao X, Bergland AO, Behrman EL, Gregory BD, Petrov DA, Schmidt PS. 2016. Global Transcriptional Profiling of Diapause and Climatic Adaptation in Drosophila melanogaster. Mol Biol Evol 33:707–720.



**Background + Objecives**



* Fruit flies occur across a latitudinal gradient, with northern populations experiencing more seasonality. 
* Diapause is a life history strategy characterized by inactivity and metabolic arrest—driven environmental cues including nutrient availability, temperature, photoperiod
* Diapause varies across latitude



==Goal is to understand the mechanisms that allow organisms to adapt to distinct environments==

* latitudinal gradients are great because environments differ continuously along this axis. 
* But! Not only can environments vary through space…but also time, especially, within a season. 
* Fruit flies originated in tropics of Africa, but have expanded world wide. And there are signatures of clinal varatiion in thse popualtions from australia, NA, SA, Africa. 
* Not only is there clinal variation in SNPS, but there are also seasonal SNPs too. 
* Northern lat SNPs are similar to those from SNPs associated with summer. (and vice versa)
* So you have SNPs resonding to the environment through space and time…..soooo

**Focus: Diapause phenotype**

* Induced by photoperiod and moderately low temperatures
* results in halt in reproduction
* increased lipid storage
* increased stress tolerance

Diapause varies along cline. (90% in temperate, 30% nootropics)



Fundamental lief history trade-off: somatic maintenance and reproduction. 



​```mermaid
graph TD
C(Benign conditions) --> I[High reproduction]
C --> T[Low survival]
E[Stressful unfavorable conditions] --> R(Higher survival)
E --> V[Low reproduction]
R --> A{Diapause}
V --> A
T--> B{Non-diapause}
I --> B
​```

Candidates for diapause: 

* Dp110
* timeless
* cough potato
* circadian clock
* insulin and insulin like growth factor



## Questions 

1. What are the common underlying mechanisms of adaptation to the environment?  (gathered from intro)
2. What are the genes and transcripts that are diff expressed as a function of diapause phenotype? 
3. Are the genes previously id as being associated with variance in diapause, diff regulated in diapause vs non diapause
4. Are the DE genes segregating for SNPs that vary in frequency with latitude and season? 
5. Is there tissue specificity? 



## Hypotheses

1. Diapause is one major determinant of adaptation to spatial and temporal environmental heterogeneity in fruit flies

Prediction:

If the natural variation of diapause does play critical roles in adaptation to environmental heterogeneity, genes differentially regulated as a function of the diapause phenotype are likely under spatially and/or temporally varying selectively pressures, thus genetic polymorphisms on these genes or in vicinity of these genes are likely to show clinal and seasonal patterns as they may have distinct cis-regulatory functions that regulates gene expression in diapausing nondiapausing individuals



Basically: Variation in gene expression associated with diapause should be under selection. One way to visualize this is how these genes vary across space(latitude) and time(season)

2. Maternal epigenetic effects: They hypothesized that higher level regulatory mechanisms such as local sharing of regulatory elements and chromatin structure may also be involved in the regulation of diapause in D. melanogaster,
3. ​

Methods



1. Collected flies from 4 popoulations (northern) east coast us: sampled in october with 50 isofemale lines each orchard
   * Bowdoin maine,
   * shoreham vt
   * harvard mass
   * middlefield ct
2. Treatment: reared at 25 C 12L :12D; 11 C 9L:15D
   * dissected, scored oocytes as :
     * D = arrest before vitellogenesis (before stage 8)
     * ND = vitellogenin was observed (stage 8 or later)
3. Measured gxp for heads and ovaries ( 2 tissue types)
   * 92 flies in D_head or O_ovary; 106 flies pooled into ND_head_ND_ovary
4. Sequencing: 100bps short reads
5. Used EBSeq (empirical baysian approach) to call genes and isoforms for gene expression. Measured RPKM ( Reads Per Kilobase of transcript per Million mapped reads)
6. Analaysis: vitellogenesis and oogenesis excluded from analysis. 
7. Location dependent expression: Looked at autocorrelation in the log2 fold changes along each chromosome (used moran's I as the metric of autocorrelation). This was compared to a randomly generated chromosome arm.
8. Test if DE genes are enriched in season/clinal sets

​```mermaid

graph TD
a(Bowdoin, Maine)--- A
b(Shoreham, VT) --- A
c(Harvard, MA) --- A
d(Middlefield, CT) --- A
A{50 isofemale lines for each population} --> B(Treatment) 
B -->|11C,9L:15D|C(Check ND or D)
C --> D[Dissect head or ovary]
D --> E[RNA-seq]
E --> G[Autocorrelation]
E --> N[Differential gene expression]
N --> V[Functional enrichment]

​```





## Results

**Fig 1.** Shows log2 fold change (D over ND) vs RPKM. So higher values of log2 fold change = diapausing individuals having more gxp. This is a plot to simply show no relationship between how highlly expressed a gene is  with expression differences. 

**thoughts**

Bad figure. What happened to the populations? That is the more interesting result.



**Table 1. Summary of gene expression level differences**

Number of genes that are differentially expressed between D and ND for head or ovary.



**Table 2. Isoform level differential expression**

Isoforms diff expressed between D and ND for head or ovary. 



**Fig2. Venn diagram of overlap between DE of genes or isoforms for head and ovary.**

* Higher number of genes with at least 1 DE isofrom in both head and ovary.
* Moderate amount of overlap
* Lowest amount of just overall DE at the gene level

**Table 3. Enriched KEGG pathways.**

* Head has more pathways than ovaries
  * more genes downregulated
* ​

**Table 4. Gene-level expression** 

Candidate genes were not differentially expressed. Why is this? 

1. cpo
2. tim
3. dp110

**Table 5. Transcript level expression**

Some of the cp, tim, and dp110 isoforms are differentially expressed. 



**Fig3. Density as a function of moran's I (autocorrelation) for each chromosome and each heard or ovary.**

* The distributions represent the null expectation in the correlation of gxp due to position on the chromosome
* It looks like all but X head chromosomes have autocorrelated gxp. 

**Table 6. Summary of overlapping genes DDE and clinal /seasonal genes**

* lots of clinal genes are differentially expressed
* less seasonal genes that are differentially expressed



**Fig 4. Enrichment of DE genes that are clinal are seasonal.**

Only downregulated in the head were enriched under clinal and seasonal. 



Take home

* Tissue specific and isoform specific expression can vary within an organism. 
* Downrgulated genes in the head are the types of genes that are related to latitude and seasonality. 
* Gene expression is autocorrelated. 

------
<div id='id-section11'/>
### Page 11: 2017-02-13. Week 5, Day 8, RNA-seq mapping and paper discussion: Johnston et al. 2016, *Molecular Ecology*

## Info update- Transcriptomics : Lauren Ash



Glossary:

1. sequence coverage: the average number of reads that align/cover known reference bases
2. read depth: total number of bases sequenced/aligned at a given reference base position
3. statistical noise: unexplained variation/randomness
4. power: probability of rejecting a false null hypothesis
5. Biological variation: natural variation in gene expression measurements due to environmental or genetic differences



1. **Background**
   * can measure differential gene expression (within population and among )
   * topics: diseases resistance, mating behavior, adaptive signfiicance 
   * connect molecular mechanisms to phenotypic/behavioral plasticity 
   * limitations: reference genomc quality, availability of genes, expense per smaple lib prep
2. **Issues**
   * under utiilization of biological replicates
   * requiring independent library preparations 
   * doesn't include pooled samples (unless pooled samples that were replicated)
     * 23/158 studies (15%) have more than 3 biological replicates
     * derive broad biological conclusions
     * prioritize sequencing depth over replication 
   * wide dynamic range can make it noisy
     1. poisson counting error (error associated with any counting experiment)
     2. non-poisson technical variance (processing, quality, different lanes, storage)
     3. biological variance (usually the largest)
3. R exercise
   * so cool, 
4. Gerenal rules of thumb
   1. more biological replicates more than increasing depth
   2. sequence more than 10 reads of depth per transcript
      * 10-20 million mapped reads per sample is sufficient
   3. Use at least 3 biological replicates per condition
   4. conduct a pilot experiment 
      1. answer 2 questions:
         * What is the and most powerful experiment that I can afford? 
         * What is the smallest fold change I can detect?
         * 7 tools in R to estimate power



## Coding: 

workflow:

1. clean adn evaluate reads (fastq)
2. Make and evaluate transcriptome assembly (fasta)
3. Map cleaned reads to transcriptome assembly (.sam files)
4. Extract data:
   1. read counts: number of reads that uniquely map to each gene
   2. get SNPs



**Transdecoder** predicts open reading frames.



## Beginning to map reads to reference transcriptome

* Using BWA to do this: Copying script

​```UNIX
cp /data/scripts/bwaaln.sh .
​```

* waht does the script look like? 



​```
#!/bin/bash 

# To run from present directory and save output: ./bwaaln.sh > output.bwaaln.txt 

myLeft='20_5-14_H_0_R1_clean_paired.fa'
echo $myLeft

myRight=${myLeft/_R1.fq.gz_left/_R2.fq.gz_right}
echo $myRight

myShort=`echo $myLeft | cut -c1-11`
echo $myShort

# bwa index /data/project_data/assembly/longest_orfs.cds  # This only needs to be done once on the reference

bwa aln /data/project_data/assembly/longest_orfs.cds /data/project_data/fastq/cleanreads/$myLeft > $myLeft".sai"
bwa aln /data/project_data/assembly/longest_orfs.cds /data/project_data/fastq/cleanreads/$myRight > $myRight".sai"
bwa sampe -r '@RG\tID:'"$myShort"'\tSM:'"$myShort"'\tPL:Illumina' \
        -P /data/project_data/assembly/longest_orfs.cds $myLeft".sai" $myRight".sai" \
        /data/project_data/fastq/cleanreads/$myLeft \
        /data/project_data/fastq/cleanreads/$myRight > $myShort"_bwaaln.sam"

​```

* Options for bwa align



​```
 bwa aln

Usage:   bwa aln [options] <prefix> <in.fq>

Options: -n NUM    max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
         -o INT    maximum number or fraction of gap opens [1]
         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]
         -i INT    do not put an indel within INT bp towards the ends [5]
         -d INT    maximum occurrences for extending a long deletion [10]
         -l INT    seed length [32]
         -k INT    maximum differences in the seed [2]
         -m INT    maximum entries in the queue [2000000]
         -t INT    number of threads [1]
         -M INT    mismatch penalty [3]
         -O INT    gap open penalty [11]
         -E INT    gap extension penalty [4]
         -R INT    stop searching when there are >INT equally best hits [30]
         -q INT    quality threshold for read trimming down to 35bp [0]
         -f FILE   file to write output to instead of stdout
         -B INT    length of barcode
         -L        log-scaled gap penalty for long deletions
         -N        non-iterative mode: search for all n-difference hits (slooow)
         -I        the input is in the Illumina 1.3+ FASTQ-like format
         -b        the input read file is in the BAM format
         -0        use single-end reads only (effective with -b)
         -1        use the 1st read in a pair (effective with -b)
         -2        use the 2nd read in a pair (effective with -b)
         -Y        filter Casava-filtered sequences

​```

* Run default parameters 
  * Only thing that would change is altering the maximum SNP differences  between mapping sequence to reference sequence
  * ​
* ==Should end up with `.sai` files!!!==

output:

​```
ls
20_5-14_H_0_bwaaln.sam              bwaaln.sh
20_5-14_H_0_R1_clean_paired.fa.sai  trim_example.sh
20_5-14_H_0_R2_clean_paired.fa.sai
​```



* ==What the header should look like:  `Transcript_Contigs::Read_ID==`

​```
head 20_5-14_H_0_bwaaln.sam 
@SQ	SN:TRINITY_DN37_c0_g1::TRINITY_DN37_c0_g1_i1::g.1::m.1	LN:303
@SQ	SN:TRINITY_DN120_c0_g2::TRINITY_DN120_c0_g2_i1::g.2::m.2	LN:381
@SQ	SN:TRINITY_DN125_c0_g1::TRINITY_DN125_c0_g1_i1::g.3::m.3	LN:642
@SQ	SN:TRINITY_DN125_c0_g2::TRINITY_DN125_c0_g2_i1::g.4::m.4	LN:528
@SQ	SN:TRINITY_DN159_c0_g1::TRINITY_DN159_c0_g1_i1::g.5::m.5	LN:696
@SQ	SN:TRINITY_DN159_c0_g1::TRINITY_DN159_c0_g1_i1::g.6::m.6	LN:396
@SQ	SN:TRINITY_DN191_c0_g1::TRINITY_DN191_c0_g1_i1::g.7::m.7	LN:309
@SQ	SN:TRINITY_DN192_c0_g1::TRINITY_DN192_c0_g1_i1::g.8::m.8	LN:318
@SQ	SN:TRINITY_DN192_c0_g2::TRINITY_DN192_c0_g2_i1::g.9::m.9	LN:321
@SQ	SN:TRINITY_DN293_c0_g1::TRINITY_DN293_c0_g1_i1::g.10::m.10	LN:315
​```

* ==What the tail should look like: `vim tail_version_sam.txt`==

​```
J00160:63:HHHT2BBXX:4:2228:14296:49089	77	*	0	0	*	CTCTGCCCCGACGGCCGGGTATAGGCGGCACGCTCAGCGCCATCCATTTTCAGGGCTAGTTGATTCGGCAGGTGAGTTGTTACACACTCCT	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJFFJJJJJJFFJJJJJJJJ7FJJJJJJJJJFJJF	RG:Z:20_5-14_H_0
J00160:63:HHHT2BBXX:4:2228:14296:49089	141	*	0	0	*	TCGGAATCCGCTAAGGAGTGTGTAACAACTCACCTGCCGAATCAACTAGCCCTGAAAATGGATGGCGCTGAGCGTGCCGCCTATACCCGGC	FJJJAJJJJJJJJFJJFFJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJ<FJJJJ-FJJJJJJJJJJAJJJJJJJFJJJJJAJJJJJF	RG:Z:20_5-14_H_0
J00160:63:HHHT2BBXX:4:2228:15026:49089	77	*	0	0	*	TTTTTCGTCACTACCTCCCCGTGTCGGGAGTGGGTAATTTGCGCGCCTGCTGCCTTCCTTGGATGTGGTAGCCGTTTCTCAAGCTCCCTC	JJJJJJJJFJJJJFJAJJJJFAJJJAJJAJ7FJJ<AJJFJFAJFAFJ<<JFAJJFJJJJJAF<AJFFJ-FJFJJFJFAJJJJFJJJFJJF	RG:Z:20_5-14_H_0
J00160:63:HHHT2BBXX:4:2228:15026:49089	141	*	0	0	*	GGTTCGATTCCGGAGAGGGAGCTTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCACTCCCGACACGGGGAGG	JJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	RG:Z:20_5-14_H_0
J00160:63:HHHT2BBXX:4:2228:22008:49089	83	TRINITY_DN30310_c1_g10::TRINITY_DN30310_c1_g10_i1::g.7248::m.7248	168	60	91M	=	70	-189	CCTCGCTCCCCGGGCGAAAGGGAATCGGGTCAATATTCCCGAACCCGGAGGCGGAGCCCTCCGTCTTCGTGGCGGTCCGAGCTGTAAAGCG	JJFJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJ	RG:Z:20_5-14_H_0	XT:A:U	NM:i:SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:91
J00160:63:HHHT2BBXX:4:2228:22008:49089	163	TRINITY_DN30310_c1_g10::TRINITY_DN30310_c1_g10_i1::g.7248::m.7248	70	60	91M	=	168	189	CCATGTGAACAGCAGTTGTACATGGGTCAGTCGATCCTAAGCCCCAGGGAAGTTCCGCTCCGAGCGGAGGCGGGCGCCCCTCTCCATGTGA	FJJJJJJJJJJJJJJFFJJJJJJJJJFJJJJJJJJJJJFFJJJJJJJJJJFJFJJJJJJFFJJJJJJJJJJJJJJJJJJJFJJJJJJJFJA	RG:Z:20_5-14_H_0	XT:A:U	NM:i:SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:91
J00160:63:HHHT2BBXX:4:2228:24647:49089	77	*	0	0	*	ACGGGCGATGTGTGCGCATTCTAGGGCTTTGAGTTGTTCATGGGCATTTTCTTTTGCTCATTACTGCTGAATCCTGTTTCAAATGGGGCTA	JFJJJJFFJJJJJJJJJFJJJJJJJJJJFJJAJJJJJJJFJJJJJJJFJJJJJJJJJF<JJJJJJJFAFJJFFFJJJJJJJJJJJJJ<JJA	RG:Z:20_5-14_H_0
J00160:63:HHHT2BBXX:4:2228:24647:49089	141	*	0	0	*	GGATAAGTGAGCTACAATCATAAATATAAGAATAAAAATATGTATGAATAATGAACTGATAGCCCCATTTGAAACAGGATTCAGCAGTAAT	AAJJJJ<JJJJFFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJFFJJJJJJJJJJ<-FFFFFAFFFAFJJJJJJJJJJJFAJFJJFJJJJ	RG:Z:20_5-14_H_0
J00160:63:HHHT2BBXX:4:2228:30168:49089	77	*	0	0	*	TCTTGAAATCTGTGGGTTTCTCGTATAGTTCAATTACAACAGGTCCTGGTTTCAACTCGTCCATTTCCATGAAGGCAAAACACTTGGTGCT	JJJFFJJJJJJJJJJJAJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJFJJJFFAFAJFJJJJFFFJJJJFFJJJJFJJJFJJ	RG:Z:20_5-14_H_0
J00160:63:HHHT2BBXX:4:2228:30168:49089	141	*	0	0	*	GAGTTTAAGCATTTCAAAGTGAAAAAGCGCACTATCAGCACCAAGTGTTTTGCCTTCATGGAAATGGACGAGTTG	JJJJJJJJJJJJJFJJAJJJJJJJJJJJJJJJJFJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFAF	RG:Z:20_5-14_H_0

​```

* Melissa code hack suggestion: 
  1. `screen`
  2. `bash bwaaln.sh`
  3. control + a + d

------
<div id='id-section12'/>
### Page 12:
------
<div id='id-section13'/>
### Page 13:
------
<div id='id-section14'/>
### Page 14:
------
<div id='id-section15'/>
### Page 15:
------
<div id='id-section16'/>
### Page 16:
------
<div id='id-section17'/>
### Page 17:
------
<div id='id-section18'/>
### Page 18:
------
<div id='id-section19'/>
### Page 19:
------
<div id='id-section20'/>
### Page 20:
------
<div id='id-section21'/>
### Page 21:
------
<div id='id-section22'/>
### Page 22:
------
<div id='id-section23'/>
### Page 23:
------
<div id='id-section24'/>
### Page 24:
------
<div id='id-section25'/>
### Page 25:
------
<div id='id-section26'/>
### Page 26:
------
<div id='id-section27'/>
### Page 27:
------
<div id='id-section28'/>
### Page 28:
------
<div id='id-section29'/>
### Page 29:
------
<div id='id-section30'/>
### Page 30:
------
<div id='id-section31'/>
### Page 31:
------
<div id='id-section32'/>
### Page 32:
------
<div id='id-section33'/>
### Page 33:
------
<div id='id-section34'/>
### Page 34:
------
<div id='id-section35'/>
### Page 35:
------
<div id='id-section36'/>
### Page 36:
------
<div id='id-section37'/>
### Page 37:
------
<div id='id-section38'/>
### Page 38:
------
<div id='id-section39'/>
### Page 39:
------
<div id='id-section40'/>
### Page 40:
------
<div id='id-section41'/>
### Page 41:
------
<div id='id-section42'/>
### Page 42:
------
<div id='id-section43'/>
### Page 43:
------
<div id='id-section44'/>
### Page 44:
------
<div id='id-section45'/>
### Page 45:
------
<div id='id-section46'/>
### Page 46:
------
<div id='id-section47'/>
### Page 47:
------
<div id='id-section48'/>
### Page 48:
------
<div id='id-section49'/>
### Page 49:
------
<div id='id-section50'/>
### Page 50:
------
<div id='id-section51'/>
### Page 51:
------
<div id='id-section52'/>
### Page 52:
------
<div id='id-section53'/>
### Page 53:
------
<div id='id-section54'/>
### Page 54:
------
<div id='id-section55'/>
### Page 55:
------
<div id='id-section56'/>
### Page 56:
------
<div id='id-section57'/>
### Page 57:
------
<div id='id-section58'/>
### Page 58:
------
<div id='id-section59'/>
### Page 59:
------
<div id='id-section60'/>
### Page 60:

------
