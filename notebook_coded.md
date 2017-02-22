01/25/2017

<u>Outline

1. What are QTN??
2. Quantitative Genetic Theory of Adaptive Traits
   - [ ] Va
   - [ ] h^2</u>
3. Methods
   - [ ] Linkage mapping
   - [ ] GWAS: "fancy regression"
   - [ ] selection scans: reverse genetic method

QTN: "Quantitative Trait Nucleotides"

- [ ] flowering time: Continuos/ Quantitative traits
- [ ] flower colors: discreet/ Mendelian traits
- [ ] thermal tolerance
- [ ] venom potency
- [ ] altitude tolerance
- [ ] defense compounds
- [ ] toxin tolerance
- [ ] draught tolerance
- [ ] altitude tolerance

Quantitative traits: you can decompose in genotypes, each genotype contributes an ammount to the phenotype. 

Fisher: Small effect vs. long effect mutations (see notebook)

# Paper discussion:

Rockman- are QTNs useful?

3 points presented: 

1. small effect are the vast majority, while the allelic frequency we are drawing from are not the most representative. "unmeasurable interesting traits"
2. ​
3. small effect alelles do not operate in the same manner: bias in effects. 



impacts of the large effect and small effects



<u>Lee et al. 2014</u>

Still useful, we might get to the small-effect QTNs with these methods

Discussion:

Use in combination with other methods

Usefulness in id. If QTN's for pop structure, and the genetic basis of an adaptation that may not apply to different populations. 

Is it worth it for the scope of a PhD program?

Communication between findings and researchers that will use this knowledge.



**01/30/2017**

##### Group projects

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



*Write up ideas into 1 page proposal using guideline provided and integrating feedback from the discussion. Include specific samples. For Monday, MS-Word and e-mail to Steve and Melissa*



Wednesday: 4 library prep types Blitz

_____________

02/01/2017

Announcements

- Send link to github online notebook to adnguyen@uvm.edu, 2 "notebooks" notes and coding
- Sign-ups
- Project proposal due next Monday by email: one per group
- Transcriptomics next week




<u>WGS</u>

- everything


<u>RNA-seq</u>: gene space expressed

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



<u>Amplicon Sequencing</u>: one gene

Def: Targeted approach for analyzing genetic variation in specific genomic regions (Amplicon= targeted gene region to be amplified via PCR w/ specific primers)

- Methods

  - Libray Prep: 

    - Extract DNA (specific region)
    - PCR 1 (16S) rRNA: Amplify gene, specific primers. —> Clean (use a gel) and can extract through the gel or use column and sequence if ok if not: —> PCR 2: barcodes and adaptors —> Clean again —> Pool —> Sequence

  - Sequencing

    - Platforms: 454, MiSeq (Ilumina)

  - Data Analysis (for 16S)

    -  Trim adaptor sequence (e.g.. AGGTTT 7bp)
    -  Align Data: use conserved areas, tell your program the seq to use 

    ​

- Applications

  - Id species
  - one gene of interest
  - Benefits: if you know the gene you want, you save time and money,

- Limitations: 

  - Computational



<u>GBS</u>: Genotyping by sequencing at the same time (before this you didn't have to develop primers and then sequence). between RNA-seq and Amplicons

"Reduced Representation"

<u>RAD seq</u>: Resticton enzyme assisted DNA sequencing

- Lots of individuals from 

- lots of SNPs across the genome

- Don't care about specific genes

  _______________

  ​


06/Feb/2017

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



______

02/08/2017

#### Info update: Transcriptomics (Lisa C.)

<u>Outline</u>

- Introduction

  Look at wild systems: 

  - non-model (no ref. genome) / non-traditional (hasn't been used traditionally) model organism
  - silent genes responding to multiple stimuli
  - novel transcripts wothout homologs in closely related model organisms

  ​

- Brief Overview of Transcriptomics

  2 Methods:

  | Microarray                            | RNA-seq                            |
  | ------------------------------------- | ---------------------------------- |
  | easy for ecological analyses          | Genome wide ecological transcripts |
  | need a reference genome (need probes) | More involved analyses             |

  Questions: 

  1. How much variation is there in gene expression & how is it structured?
  2. How do environmental stimuli affect gene expression?
  3. How does gene expressiona affect phenotype?

- Main Questions

  1) Evolutionary processes

  - Gene expression heretable— natural selection
  - Qst - Fst comparissons
  - eQTL expression quantitative trait loci mapping, expression of populations within populations
  - Epigenetics
  - Macroevolution: drift, selection, bottleneck

  2) Q2: environmetal stimuli

  - abiotic stress
  - environmental heterogeneity
  - host-parasite interactions
  - selective biotic and abiotic interactions: at levels: molecular, genotype, phenotypic plasticity
  - a snapshot at that time of expression: flash freeze to prevent environmetal effects (if not RNA later)
  - time course analysis: transcritional response through time ($$)

  3) Q3: Gene expression and phenotype

  - Alternate phenotypes 
    - moving from correlation to causation: transgenics, RNAi, CRISPR/CAS (use of knockdown of genes for this analysis)

- Future Directions

  - Use of combined Microarrays and RNA-seq to test (need availabe probes depending on the questions)
  - Database for proposed ecological annotations
  - Address problems
    - Bias in signal: in commonly expressed genes
    - heterologous arrays: just microarrays
    - Polyploidy
    - RNA Pooling
    - Statistical Analysis
    - Unannotated Genes

Glossary:

Transcritomics: the srudy of the transcriptome which is the complete set of RNA transcripts produces by the genome

Qst- Fst comparissons: A means to distinguish between genetic drift and natural selection in driving differentiation in a population. (Qst= amount if variation in quantitative traits in a popultaion) (Fst= variation in a neutral loci)



_______

02/13/2017

##### <u>RNA-seq Info Update 3</u>

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



##### <u>Paper discussion</u>

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
- ​


_______

02/15/2017



Info Update: SNPs & Population Genomics



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

##### Paper discussion : led by Andrew

clinal variation: varies by a gradient, in this case 

diapause: developmental quiescence, inactive stage. could happen on any life stage:

​	caused by: temperature variation, photoperiods

Are these SNPS associated with this diapause?

Methods:

50 isofemale lines (maintain genotype)

Treatment: decreased temp and photoperiod

Replicates: Diapause, NonDiapause 2 tissues: head and ovary

isoform: splice variant of the same gene

Moran's I = how samples that are close in space are functionally similar. 

Results:

Fig 1: dif in gene expression between D and ND between each tissue. 

RPKM = FPKM: normalized read by lenght of the transcript, and how much you sequenced was sampled. scale the reads by two factors

Fig 3: DE genes (patterns) associated with their location in chromosome neighborhoods. Maybe related to chromatin structure?

enrichment: in which genes is it enriched in...





