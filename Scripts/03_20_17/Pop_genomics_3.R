# 03_20_17 EcoGenomics LAB

# Set your working directory to where you downloaded your results files:


# List the files in this directory -- you should see your results output from VCFTools if the download was successful
#shows all files in the working directory
list.files()

# Let's do the allele freq comparisons first:
H_freq <- read.table("H_AlleleFreqs.frq", header=T)
S_freq <- read.table("S_AlleleFreqs.frq", header=T)

# Since these files have identical numbers of SNPs in the exact same order, we can concatenate them together into one large dataframe:
All_freq <- merge(H_freq, S_freq, by=c("CHROM", "POS"))

# Check the results of your merge to make sure things look OK
str(All_freq) # shows the structure of the data, make sure you edited the headers
head(All_freq)

# Looks good, now let's calculate the difference in minor allele frequency at each SNP and plot as a histogram
All_freq$diff <- (All_freq$H_ALT - All_freq$S_ALT)

hist(All_freq$diff, breaks=50, col="red", main="Allele frequency difference (H-S)")

# Looks like most loci show little difference (i.e., likely drift), but perhaps a few show very large differences
# between healthy and sick (drift or selection?)

# How do these highly divergent frequenices compare to Fst at the same SNPs?
fst <- read.table("HvS_Fst.weir.fst", header=T)
str(All_freq.fst)
All_freq.fst <- merge(All_freq, fst, by=c("CHROM", "POS"))

plot(All_freq.fst$diff, All_freq.fst$WEIR_AND_COCKERHAM_FST, xlab="Allele frequency difference (H-S)", ylab="Fst", main="Healthy vs. Sick SNP divergence")

# Which are the genes that are showing the highest divergence between Healthy and Sick?
# within ALL_freq.fst data frame, ask which have the Fst>0.2
All_freq.fst[which(All_freq.fst$WEIR_AND_COCKERHAM_FST>0.2),]

# to know the Fst for a particular transcript
# Asked all rows for this transcript ID
All_freq.fst[which(All_freq.fst$CHROM=="TRINITY_DN42225_c1_g1_TRINITY_DN42225_c1_g1_i1_g.12458_m.12458"),]

# Add a point to the plot. first I need to know the x and y coordinate to a point
# points then the axis (x,y) for the SNP
points(0.2500000, 0.2415350, col="purple", cex=4)
