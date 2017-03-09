# 2017 Ecological Genomics Course

### Author: Laura Caicedo-Quiroga

## Overall Description of notebook

Here you can see my online notebook for the Ecological Genomics course offered for the Spring of 2017. I will take notes on themes on the subject, paper discussions and methods discussed on the lab portion of the class. 

## Date started: (2017-01-18)

## Date end:   -



<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.  

### Table of contents

- [Page 1: 2017-01-18](#id-section1). Introduction to Ecological Genomics
- [Page 2: 2017-01-23](#id-section2). Day 2 (Info update rubric, course notes, paper discussion)
- [Page 3: 2017-01-25](#id-section3). Day 3
- [Page 4: 2017-01-30](#id-section4) . Day 4: Group projects
- [Page 5: 2017-02-01](#id-section5). Day 5
- [Page 6: 2017-02-06](#id-section6). Day 6
- [Page 7: 2017-02-08](#id-section7). Day 7
- [Page 8: 2017-02-13](#id-section8). Day 8
- [Page 9: 2017-02-15](#id-section9). Day 9
- [Page 10: 2017-02-22](#id-section10). Day 10
- [Page 11: 2017-02-27](#id-section11). Day 11
- [Page 12: 2017-03-01](#id-section12). Day 12
- [Page 13: 2017-03-06](#id-section13). Day 13
- [Page 14: 2017-03-08](#id-section14). Day 14
- [Page 15 2017-03-08:](#id-section15). Homework #2
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

### Page 1: 2016-01-18. Day 1, Ecological genomics: Introduction

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

### Page 4: 2017-01-30. Day 4

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

### Page 5: 2017-02-01. Day 5

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

### Page 6: 2017-02-06. Day 6

#### <u>Info update: RNA-seq</u> 

Prior Set-up:

1. Approach
2. Experimental design
3. Library prep
4. Sequencing (Facility)
5. Receive data
6. Computer/ Server setup

<u>RNA-seq Workflow</u>

1) Clean reads [.fastq, .fq, .gz]

- adaptors
- nucleotide quality
- length

2) Evaluate quality

3) de novo transcriptome assembly (making reference) [.fasta]. If there is a reference, skip this step 

- evaluate assembly: compare to closely related species, CEG (core eukaryotic gene), N50, # of contigs
- Annotation (BLAST search)
  - BLASTx: nr (gene annotation), uniprot database (curated database of proteins): using gene name Gene Ontology (GO) functional contigs. Amino acid sequences used 

4) Map reads to reference transcriptome, alignment files generated. How many reads correspond and missalignments, decide how many to allow. Computationally intensive and important step. [.sam = sequence alignment file]

5A) Extract read count info: # of reads that map to each contain for each sample

5B) Identify SNPs

6A) DGE= Differential gene expression analyses, co-expression network

6B) Population Genomics: genetic differentiation, population structure, demographic history, test for signatures of selection.

**Discussions: Brief summary of objectives and methods**

<u>Dunning_et al._2014: Discussion</u>

New Zealand stick insects 4 species, dif population with dif. adaptations to cold environment. Objectives: genes expressed differentially after cold shock by species and populations ans species. Make predictions about climate change response, whether they share this response or not.

Methods: 

- Cold-shock, RNA extraction, 50 short paired-end reads

- Cleaning:

  - Triming of low quality bases   (Phred <30)

- Quality 

- Assembly: Transcriptome de novo assembly using Trinity 

  - seq Mapped back, unigenes identified 

- Differential Expression

  - 3 programs: DESeq, edgeR, baySeq; deal with the variance generated differently. "Genes were considered as differentially expressed if they were significant for at least one of the differential expression analyses.""

  - Enrichment analysis: gene ontology analysis with BLAST2GO

    - ef. 60 K unigenes: ~10% for cuticular HC-nr, 2K DGE ~20% cuticular HC

      |      | total | chc   |      |
      | ---- | ----- | ----- | ---- |
      | all  | 60K   | 600   | 10%  |
      | DEG  | 2K    | *400* | 20%  |


- Biological replicate (by individuals) and technical replicates (1 individual)
- Fig 2: mtDNA (quickly evolving) only not as nuclear phylogeny (28S=slow evolving). Possible introgressions within species. 
- Fig 3: Use of SNP data a) want highest resolution, lowest cross-validation error with higher level of K, K=3 selected. bars with more than one color represent admixture. b) PCA, next method to test this. Much better
- Fig 4: Sewell Peak has high quantity of up and down regulated genes for cold-response, very few shared responses. CS + CO = DR , CS and CO are muttually exclusive
- Table 2: GO
- Fig 5: 6 genes selected, changes in gene expression between control and col-shocked, red to green= down to up regulation. qPCR to validate read counts from RNA-seq

THM (Take home message):

Functional enrichment of cuticular structure genes, mostly in higher altitudes. 

_______

<div id='id-section7'/>

### Page 7: 2017-02-08. Day 7

#### <u>Info update: Transcriptomics (Lisa C)</u>

<u>Outline</u>

1. Introduction

   Look at wild systems: 

   - non-model (no ref. genome) / non-traditional (hasn't been used traditionally) model organism
   - silent genes responding to multiple stimuli
   - novel transcripts wothout homologs in closely related model organisms

   ​

2. Brief Overview of Transcriptomics

   2 Methods:

   | Microarray                            | RNA-seq                            |
   | ------------------------------------- | ---------------------------------- |
   | easy for ecological analyses          | Genome wide ecological transcripts |
   | need a reference genome (need probes) | More involved analyses             |

   Questions: 

   1. How much variation is there in gene expression & how is it structured?
   2. How do environmental stimuli affect gene expression?
   3. How does gene expressiona affect phenotype?

3. Main Questions

   1) Evolutionary processes

   - Gene expression heritable— natural selection
   - Qst - Fst comparissons
   - eQTL expression quantitative trait loci mapping, expression of populations within populations
   - Epigenetics
   - Wild populations
   - Macroevolution: drift, selection, bottlenecks

   2) Q2: environmetal stimuli

   - abiotic stress
   - environmental heterogeneity
   - host-parasite interactions
   - selective biotic and abiotic interactions
   - at levels: molecular
     - among genotypes
     - phenotypic plasticity
   - a snapshot at that time of expression: flash freeze to prevent environmetal effects (if not RNA later)
   - time course analysis: transcritional response through time ($$)

   3) Q3: Gene expression and phenotype

   - Alternate phenotypes 
     - moving from correlation to causation: transgenics, RNAi, CRISPR/CAS (use of knockdown of genes for this analysis)

4. Future Directions

   - Use of combined Microarrays and RNA-seq to test (need availabe probes depending on the questions)
   - Database for proposed ecological annotations
   - Address problems
     - Bias in signal: in commonly expressed genes
     - ~~heterologous arrays~~: just microarrays
     - Polyploidy
     - RNA Pooling
     - Statistical Analysis
     - Unannotated Genes

Glossary:

Transcritomics: the study of the transcriptome which is the complete set of RNA transcripts produces by the genome

Qst- Fst comparissons: A means to distinguish between genetic drift and natural selection in driving differentiation in a population. (Qst= amount if variation in quantitative traits in a popultaion) (Fst= variation in a neutral loci). 

#### Paper discussion by Laura Caicedo: ___

Dessert drosophila (S.America) 2 cactus hosts differ in alkaloid content and composition.

_______

<div id='id-section8'/>

### Page 8: 2017-02-13. Day 8

#### <u>RNA-seq Info Update 3</u>

Outline:

- Background
- Issues
- R exercise
- General Rules of Thumb

A) Background

- Enables DE examination (inter-population individual)
  - disease resistance
  - mating behavior
  - adaptive significance
- Molecular mechanisms  
  - phenotypic/behavioral plasticity
  - migration patterns


- Limitations
  - reference genome quality
  - gene annotation availability
  - expense per sample library prep.

B) Issues

- Under utilization of biological replicates
  - Requiring independent library preparations
  - doesn't include pooled samples
  - 23/258 studies (15%) > 3 biological replicates
  - derive broad biological conclusions
- Prioritize sequence depth over replication (one of the key trade-offs), **big problem**
- Wide dynamic range of RNA-seq data —> Noisy
  - Poisson counting error
  - non-poisson technical variance (associated with storage of sample...)
  - Biological variance (natural), ussually the largest variance

C) R exercise

- Effect size measured in fold change: 
  - control gene: 15 —> 30, fold change= 30/15=2 (up regulation)
- lower variance gives higher power
- with high biological variance, you don't know where the variance is coming from therefore power is low. 

D) General Rules of Thumb

- Use more biological replicates to estimate biological variance more effectively
- Sequence depth > 10 reads/transcript (sufficient enough to minimize technical error).  ~ 10-20 million mapped reads/ samples
- 3 Biological replicates per condition
- conduct a pilot experiment (to determine number of replicates)
  - What is the best/powerful experiment I can afford?
  - What is the smallest foldchange I can detect?

<u>Glossary</u>

1. Sequence coverage = the average # of reads that align/ "cover" known reference bases
2. Read depth = total # of bases sequenced/aligned at given reference base position.(multiple reads covering the exact same base pair in order to have depth). Low expression genes, may not be detected with low read depth. 
3. Statistical noise = unexplained variation/ randomness 
4. Power = probability of refecting false null hypothesis
5. Biological variation = natural variation in gene expression measurements due to environmental or genetic diffs. 

#### <u>Paper discussion</u>

Johnston et al. 2016: Seasonal gene expression in a migratory songbird

Measure of transcriptome wide gene expression on migratory Swanson's thrushes (2 subspecies). C. u. ustulatus- coastal and C. u. swainsoni- inland



Methods: 

- Collected in summer using mist nets. sampling (sacrifice) on 2 seasons (summer and autumn simmulated). Migratory and non-migratory states: nocturnal activity for migratory birds (all post reproductive). Biological replicates all higher than 3. 

- gene expression: still want to use unpaired reads, weigh them equally as paired-end reads even if they don't match reference genome 

- q-value as opposed to P-value, they have too many genes and regressions

  ​

Results:

- Fig 1: A= "scatterplot" y axix is Q-values (predicted # of false positives hanging on the tail, p-value distribution and find false discovery rate) here when higher to 0.2 and closer to 0 it is significant. x-axis: log2 (foldchange of 2). red for genes highlighted in B. B=
- Fig 2

THM:

- expected a set of genes fo be DE for circannual changes but found others, related to light changes
- Genes as hubs of gene expression networks

______

<div id='id-section9'/>

### Page 9: 2017-02-15. Day 9

#### <u>Info Update:</u> SNPs & Population Genomics

SNP data- from expressed sequences: To overcome drawbacks of using genomic data ( is this functionally relevant?)

<u>Methods</u>

1) Tissue: breadth of tissue from multiple development stages (captures variation). Issue = exon skipping.

2) Pool & Sequence libraries: ~ 30-100m paired-end long reads

3) Process raw sequence data: important for SNP detection

4) Digital Normalization: to remove high coverage reads & associated errors. Reduces sampling variation, and no quantitative info. 

5) Assemble cleaned paired long reads 

6) Prune: reduce DNA contamination, non-coding RNA, gen fragments

7) Assembly evaluation: If Reference Genome use, if not use COGS

8) Next: SNP detection

Software that looks for constant patterns of sequence variation (to reduce errors without removing SNPS)

Problems:

- Sequence error : eliminated SNPs of low frequency
- Artifacts caused InDels: filter SNP clusters near INDELS, can also use Quality scores

9) Validation:

Primers, sequencing, mass spectrometry

<u>Applications:</u>

- Different in population structure
- Natural Selection acting on a particular loci

<u>Outliers:</u> for a given locus, what is the level of differentiation compared to difference across genome? (use outliers)

- on a graph of frequency vs. divergence Fst: a normal distribution, on tail there will be directional selection. (really high values of Fst)

<u>Non-outliers</u>

- test high Fst loci for other features associated with selection
  - fitness 
  - enrichment for certain funtional roles

— What to do if you don't have mult. populations? use other footprints (not Fst) 

<u>Glossary</u>

SNPs: Single nucleotide polymorphism single base difers between 2 genomes

InDels: Insertion/deletions, single base has been deleted/inserted into one genome relative to another.

Fst: % of genetic material explained by diggerences among populations (in this case all loci). only if you have mult pops. 

#### Paper discussion : led by Andrew

- clinal variation: varies by a gradient, in this case 
- diapause: developmental quiescence, inactive stage. could happen on any life stage:
- ​caused by: temperature variation, photoperiods
- Are these SNPS associated with this diapause?

<u>Methods:</u>

1. 50 isofemale lines (maintain genotype)
2. Treatment: decreased temp and photoperiod
3. Replicates: Diapause, NonDiapause 2 tissues: head and ovary
4. isoform: splice variant of the same gene
5. Moran's I = how samples that are close in space are functionally similar. 
6. Results:

Fig 1: dif in gene expression between D and ND between each tissue. 

RPKM = FPKM: normalized read by lenght of the transcript, and how much you sequenced was sampled. scale the reads by two factors

Fig 3: DE genes (patterns) associated with their location in chromosome neighborhoods. Maybe related to chromatin structure?

enrichment: in which genes is it enriched in...

______

<div id='id-section10'/>

### Page 10: 2017-02-22. Day 10

#### Paper Discussion: Dixon et al. 2017

Fig1: experimental design and results. E &F: Dam's effect is greater than sire's

2 hypotheses: maternal efffect or biparental expression effect (by epigenetic modification of the nuclear genome)

Fig2: Gene expression before heat treatment. resistance associated genes vs. tolerance associated genes. C: dendrogram= shows relatedness of genes. 

Fig3: tolerance vs. stress: opposite expression

___

<div id='id-section11'/>

### Page 11: 2017-02-27. Day 11

**  I was absent this day **

Guest speaker: DR. Scott Edwards

Paper Discussion: Reticulation, divergence, and the phylogeography–

phylogenetics continuum - Edwards et al.

_____

<div id='id-section12'/>

### Page 12: 2017-03-01. Day 12

#### Plan

- DESeq2 wrap up (lab)
- Homework assignment: H/S on Int. vs Sub. and another ignoring where they came from. (for entire data)
  - Verbal description: nuts and bolts of what we have done, include models
  - part of electronic notebook dedicated to assignment and refer to it on hw. 
  - from seq data to expression counts
  - due: next wednesday by end of day, word doc, email Steve, Melissa and Andrew
- Info. Update

#### <u>Info Update:</u> WGCNA

<u>Outline</u>

1. Overview of WGCNA
2. Network construction
3. Module Detection
4. Incorporation of external info
5. Topological properties
6. Other features
7. Limitations

— 

1. WGCNA = Weighted gene correlation network analysis
   - R package: apply correlation network methods to describe correlation (co-expression) patterns among genes in micro-array samples
2. Network construction -> Module identification -> Relationship of modules to external information -> Relationship between/within modules -> Finding key drivers in modules of interest.

- Network construction: use of nodes (genes) correlated by their levels of gene expression, edge (strength of correlation in gene expression).

- Package provides different co-expression measures.

  - Signed networks: positive correlation in expression
  - Unsigned networks: absolute value of correlation in expression

- Matrix example: 

  ​			*i*l….m (individuals)

  ​			  1  2  3 …………………..

  l.*..n(data)	[  (  sample's genes  )  ]

  genes/traits	[ g1 g2 g3 ………………...]

- Unweighted Network analysis

  - hard threshold for wether genes are linked or network for expression
  - can loose information

- Weighted Network Analysis

  - package allows soft or hard threshold 

1. Module Detection
   - removing weak connection (like collapsing branches on a tree)

- module: unsupervised clustering of nodes (no a priori defined gene sets)
  - summarize the profiles of the modules -> involves eigengene

1. Incorporation of external info into network
   - depends on the questions
   - Gene significance: assigns a positive number with each gene

— 7. Limitations: limited to undirected networks

— Glossary:

- Eigengene: when a sample trait is available, one can correlate the module eigengenes with its outcome.

______

<div id='id-section13'/>

### Page 13: 2017-03-06. Day 13

#### <u>Info update:</u> Population Genomics (Steve)

Population Genomics: (def) Pop gen at more precise scale and ask how is selection shaping the genome.

- SNPs (many)
- sampling unit is individuals within species
- processes:
  - Population structure
  - diversity within populations
  - selection
    - positive
    - negative

<u>Pipeline:</u>

Raw reads —> clean —> assemble "draft transcriptome" —> mapped reads —> TRANSCRIPTOMICS: Count # reads/ transcripts —> DGE!

or

…mapped reads—> POPULATION GENOMICS: Call SNPs + genotypes (homozygote, heterozygote…) —> allele frequencies, SFS, (pi)

<u>Challenges of SNP calling:</u>

1) Sequencing Error (Illumina 1:100)

- Criteria for knowng if SNP is true: apply filters

  - Minor Allele Frequency (how many individuals have the SNP, discard those rare across individuals)

  - Depth (less reads covering a position), small sample size

  - Heterozygotes: true or not, therefore how to assign a probality to this

    - if G=AT, predict A=T=0.5
    - user defined threshold for probability and filter 
    - or use them (bayesian statistics), incorporate probabilities to _

  - Paralogy: gene duplicates, causes false heterozygote (all individuals in the population) you violate HWE= 1= p^2 +2pq + q^2

    —A———— and —B———— (B has a SNP *T)

Diversity: 

- Pi = pairwise nucleotide diversity. Expected heterozygosity

  - Sequences *i + j* 

    Pi= Sum  Xi Xj Piij
    Pi(syn)= 4Ne u (miu at synonimous sites)

    Pi(nonsyn)= (almost always deleterious)
    Pi(syn)/Pi(nonsyn) : ratio

Glossary:

paralog= gene duplicate

Pi= pairwise nucleotide diversity

— 

Paper discussion: Alex

Gayral et al. 2013

Expectection: small Ne will have low Pi(syn) and Pi(nonsyn) and higher ratio than larger Ne populations. (mutations are occurring at meiosis)

Paralog filtering: eliminate false heterozygotes, also you don't have a reference genome to detect where there may be paralogs. 

Methods (p.11): Pileline using Ilumina 454 (old)

- Why outgroup? Higher differences between ingroup and outgroup and within ingroup. ratios of syn and nonsyn diversity when comparing ingroup to outgroup. How many mutations have happened? (drift vs selection in the lineage)
- If no ORF: less confidence in the alignments 
- Ingroup: clean paralogs(within species) gene duplicates, outgroup: find orthologs (between species-same gene)

Fig 2: more rare things would 

rare alleles most likely in heterozygotes than homozygotes

-----

<div id='id-section14'/>

### Page 14: 2017-03-08. Day 14, Population Genomics 2

Ne : Parameter that affects everything, and everything seems to affect it

#### Info Update:  by Kattia

<u>Rate of evolution due to relationships of effective pop size and substitution rate (NeRR)</u>

Conventions:

- Nc= census pop size
- u = mutation rate
- s = selection coefficient
- w = fitness
- Ne RR= relationship b/w NE and substitution rate
- DEF= Distribution of fitness effects (each distribution has different effects)

Outline:

1. Ne
2. Mutations
   1. overview
   2. substitutions
   3. 5 types
   4. variation
   5. u
3. NeRR + Linkage

--

1. Ne: 4 methods measure:

- from sp life history -> Ne=4NmNf/(Nm + Nf)
- from variance in allele frequency b/w genes
- from genetic polymorphism data
- correlated trait qith body size

Ne: varies across species and across the genome

- genome: due to
  - genetic hitchhiking/selective sweeps: negative mut increase because they are attached to positive (when recombination doesn't happen)
  - background selection: positive mut increase because they are attached to negative  (when recombination doesn't happen)
  - fewer sex chromosomes than autosomes)

2. Mutation can occcur:
   - whole gene or chromosome: duplication, inversion, deletion, translocation
   - base level:  (our focus) point mutation—> substitutions:
     - transitions: purine to purine, pyrimidine to pyrimidine
     - transversion: purine to pyrimidine and back
   - 2 types: 
     - synonymous: silent site —> DNA seq change= AA doesn't change: ~~NS~~ but by drift
     - non synonymous: replacement mutations—> DNA sequence change = change in AA: due to NS by purifying selection and positive selection. 
   - 5 classes: 
     - Neutral: w (effect on fitness) ~0, <1/Ne, No effect on NS, drift
       - most commonly used to estimate Ne
     - Slightly deleterious and Slightly advantageous: small effect on w, n: 1/Ne, 2n: 1/2Ne. we see effect on NS and drift
     - Deleterious: big effect on w, > 1/Ne and effect on NS, - NeRR
     - Advantageous: big effect on w, > 1/Ne and effect on NS, + NeRR
   - Variation due to
     - generation time: short gen. time-> high u: + NeRR
     - selection: lowers u
3. NeRR + Linkage
   - NeRR affected by: 
     - selective sweeps
     - clonal interference= 2 or more adaptive mutations that originate in each and compete for nex generation
4. NeRR + Spatiotemporal variation
5. NeRR + All mutations DEF
6. w landscapes: alternative app. to study NeRR: trait vs. w (curve w/ optimal peak where pops further from peak have faster muntation)

THM: The study of NeRR helps to understand the process that drives and limits evolution. Drift and selection are the most important forces determining NeRR. We need to work on better way to estimate Ne. With time, less $ m DNA seq will permit estimation of Ne, substitution rates, mutation rates. Most mutations are deleterious (Neutral theory: neutral mut are mostly deleterious).

Haploids: (theta) Ne u ~ Pi syn

Diploids: (theta) 2Ne u

Diploids (separate sex): - 4Ne u = Pi syn —> Ne= Pi/ 4 u

--

Paper discussion: Rominguier et al. 2014

Pi syn diversity across metazoa 

-------

<div id='id-section15'/>

Page 15: 2017-03-08. Homework #2

```R
# Laura Caicedo-Quiroga
# Script for HW#2

library("DESeq2")

library("ggplot2")

countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

#subset intertidal and subtidal
colDataINT <-subset(colData, colData$location=="int")
colDataSUB <-subset(colData, colData$location=="sub")

countDataINT <-countData[, which(colnames(countData) %in% row.names(colDataINT))]
countDataSUB <-countData[, -which(colnames(countData) %in% row.names(colDataINT))]
dim(countDataINT)
dim(countDataSUB)

#################### MODEL NUMBER 1: TEST EFFECT OF HEALTH 
#First for Intertidal

ddsINT <- DESeqDataSetFromMatrix(countData = countDataINT, colData = colDataINT, design = ~ health)

dim(ddsINT)
ddsINT <- ddsINT[ rowSums(counts(ddsINT)) > 100, ]
dim(ddsINT)

colData(ddsINT)$health <- factor(colData(ddsINT)$health, levels=c("H","S")) 
#sets that "healthy is the reference

ddsINT <- DESeq(ddsINT) 

resINT <- results(ddsINT)
resINT <- resINT[order(resINT$padj),]
head(resINT)
summary(resINT)

#Intertidal
plotMA(resINT, main="DESeq2INT", ylim=c(-2,2))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(ddsINT, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","day")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,500)
p

## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(ddsINT, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p
p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() + ggtitle("Intertidal")
p

#Model 2 for Subtidal

ddsSUB <- DESeqDataSetFromMatrix(countData = countDataSUB, colData = colDataSUB, design = ~ health)

dim(ddsSUB)
ddsSUB <- ddsSUB[ rowSums(counts(ddsSUB)) > 100, ]
dim(ddsSUB)

colData(ddsSUB)$health <- factor(colData(ddsSUB)$health, levels=c("H","S")) 
#sets that "healthy is the reference

ddsSUB <- DESeq(ddsSUB) 

resSUB <- results(ddsSUB)
resSUB <- resSUB[order(resSUB$padj),]
head(resSUB)
summary(resSUB)

#Subtidal
plotMA(resSUB, main="DESeq2SUB", ylim=c(-2,2))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(ddsSUB, gene="TRINITY_DN42073_c0_g1_TRINITY_DN42073_c0_g1_i1_g.12173_m.12173", intgroup=(c("health","day")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,500)
p

## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(ddsSUB, gene="TRINITY_DN42073_c0_g1_TRINITY_DN42073_c0_g1_i1_g.12173_m.12173", intgroup=(c("health","score")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p
p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() + ggtitle("Subtidal")
p

### Model 3 accounting for location
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)
dim(dds)
dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)

colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) # (factor)sets that "healthy is the reference

dds <- DESeq(dds) 

res <- results(dds)
res <- res[order(res$padj),]
head(res)
summary(res)

# ALl accounting for location
# MA-Plot
plotMA(res, main="DESeq2ALL", ylim=c(-3,3.5))


## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(dds, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,500)
p
## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(dds, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p
p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() + ggtitle("All")
p



############## PCA plots
# For all
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("score"))
plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))


# For Intertidal
vsdINT <- varianceStabilizingTransformation(ddsINT, blind=FALSE)

plotPCA(vsdINT, intgroup=c("score"))
pca2<-plotPCA(vsdINT, intgroup=c("health"))
pca2 + ggtitle("Intertidal")
plotPCA(vsdINT, intgroup=c("day"))



# For Subtidal
vsdSUB <- varianceStabilizingTransformation(ddsSUB, blind=FALSE)

plotPCA(vsdSUB, intgroup=c("score"))
plotPCA(vsdSUB, intgroup=c("health"))
plotPCA(vsdSUB, intgroup=c("day"))

```



