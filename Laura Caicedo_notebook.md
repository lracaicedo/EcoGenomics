# 2017 Ecological Genomics Course

### Author: Laura Caicedo-Quiroga

## Overall Description of notebook

Here you can see my online notebook for the Ecological Genomics course offered for the Spring of 2017. I will take notes on themes on the subject, paper discussions and methods discussed on the lab portion of the class. 

## Date started: (2017-01-18)

## Date end:   -



<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.  

### Table of contents

- [Page 1: 2017-01-18](#id-section1). Introduction to Ecological Genomics
- [Page 2: 2017-01-20](#id-section2). Day 2 (Info update rubric, course notes, paper discussion)
- [Page 3: 2017-01-25](#id-section3). Day 3
- [Page 4: 2017-01-30](#id-section4) . Day 4: Group projects
- [Page 5: 2017-02-01](#id-section5). 
- [Page 6: 2017-02-01](#id-section6). 
- [Page 7: 2017-02-03](#id-section7). 
- [Page 8: 2017-02-06](#id-section8). 
- [Page 9: 2017-02-08](#id-section9). 
- [Page 10: 2017-02-10](#id-section10). 
- [Page 11: 2017-02-13](#id-section11). 
- [Page 12:](#id-section12).
- [Page 13:](#id-section13).
- [Page 14:](#id-section14).
- [Page 15:](#id-section15).
- [Page 16:](#id-section16).
- [Page 17:](#id-section17).
- [Page 18:](#id-section18).
- [Page 19:](#id-section19).
- [Page 20:](#id-section20).
- [Page 21:](#id-section21).
- [Page 22:](#id-section22).
- [Page 23:](#id-section23).
- [Page 24:](#id-section24).
- [Page 25:](#id-section25).
- [Page 26:](#id-section26).
- [Page 27:](#id-section27).
- [Page 28:](#id-section28).
- [Page 29:](#id-section29).
- [Page 30:](#id-section30).
- [Page 31:](#id-section31).
- [Page 32:](#id-section32).
- [Page 33:](#id-section33).
- [Page 34:](#id-section34).
- [Page 35:](#id-section35).
- [Page 36:](#id-section36).
- [Page 37:](#id-section37).
- [Page 38:](#id-section38).
- [Page 39:](#id-section39).
- [Page 40:](#id-section40).
- [Page 41:](#id-section41).
- [Page 42:](#id-section42).
- [Page 43:](#id-section43).
- [Page 44:](#id-section44).
- [Page 45:](#id-section45).
- [Page 46:](#id-section46).
- [Page 47:](#id-section47).
- [Page 48:](#id-section48).
- [Page 49:](#id-section49).
- [Page 50:](#id-section50).
- [Page 51:](#id-section51).
- [Page 52:](#id-section52).
- [Page 53:](#id-section53).
- [Page 54:](#id-section54).
- [Page 55:](#id-section55).
- [Page 56:](#id-section56).
- [Page 57:](#id-section57).
- [Page 58:](#id-section58).
- [Page 59:](#id-section59).
- [Page 60:](#id-section60).

------

<div id='id-section1'/>

### Page 1: 2016-01-18. Ecological genomics: Introduction

#### **Uses**

- Understand the genetic source of adaptations
- Approaches and tools for studying relationship: <u>genome/phenotype/environment</u>
  - Ecological issues such as nutrient cycling, population structure, life history vairation , trophic interaction, stress responess, and adpatation to environmental change 
- Environmental change and adaptation studies


- On any system, not only model organisms 
- Next-gen sequencing: increasingly useful but with huge datasets (computational challenge)    

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

- Speciation, hybridization
- Local adaptations        

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

- Flow cell ,has lanes and oligos (eg. P7, P5) that attaches to your DNA that has adaptors thatsample by ligation. Your DNA is attached to a sequence adaptor and may have anID barcode. Then **bridge amplification**where the pieces bend over and amplify (to copy and make signal stronger). Bridgereleased and then cluster generation (clusters of the same sequence). For eachcluster a snapshot is taken of the nucleotides, as the polymerase is adding thenucleotides. 
- PacBio: Single Molecule in Real Time. Givesreally long reads but accuracy is lower. 


- Drawings on notebook*

6) Learning Activity

### Glossary:

Ilumina reads: short (50bp), long (100,150,300 bp), 

SMRT: extra long (10000 – 60000 bp). Can also be single vs.paired end. 

------

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

1. Quantitative Genetic Theory of Adaptive Traits
   - Va
   - h^2
2. Methods
   - QTL - Linkage mapping (forward genetics)
   - GWAS- Genome-wide association studies: "fancy regression" (forward genetics)
   - Selection scans (reverse genetics)

#### Paper discussion:

<u>Rockman, (directed by Allison Brody)</u>

- are QTNs useful?

3 points presented: 

1. small effect are the vast majority (mendelian), while the allelic frequency we are drawing from are not the most representative. "unmeasurable interesting traits"
2. ​Theory
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