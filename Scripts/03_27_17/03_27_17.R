# Set your working directory to where you downloaded your results files:


list.files() # Do you see your downloaded files there? If not, double check to make sure you've set your working directory to the right spot

# We'll need to install 2 packages to work with the SNP data:
install.packages("vcfR") # reads in vcf files and proides tools for file conversion 
install.packages("adegenet") # pop-genetics package with some handy routines, including PCA and other multivariate methods (DAPC)

# ...and load the libraries
library(vcfR)
library(adegenet)

#Read the vcf SNP data into R
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")

# The adegenet package uses a highly efficient way of storing large SNP datasets in R called a "genlight" object. The following function creates a genlight object from your vcf:
gl1 <- vcfR2genlight(vcf1)
print(gl1) # Looks good! Right # of SNPs and individuals!

# For info, try:
gl1$ind.names #names in the order of vcf files
gl1$loc.names[1:10] #first 10 locus names
gl1$chromosome[1:3] #just transcript IDs

# Notice there's nothing in the field that says "pop"? Let's fix that...
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort by Individual ID, just like the VCF file
str(ssw_meta)

# Confirm the ID's are ordered the same in gl1 and ssw_meta:
gl1$ind.names #see names and sorted from small to big
ssw_meta$Individual

gl1$pop <- ssw_meta$Location # assign locality info
gl1$other <- as.list(ssw_meta$Trajectory) # assign disease status, into the "other" slot

# WE can explore the structure of our SNP data using the glPlot function, which gives us a sample x SNP view of the VCF file
glPlot(gl1, posi="bottomleft") #generates a "heatmap" for every SNP where I have data, and where it is missing (white)

# Now, let's compute the PCA on the SNP genotypes and plot it:
pca1 <- glPca(gl1, nf=4) # nf = number of factors: of PC axes to retain (here, 4)
pca1 # prints summary
pca1

# Plot the individuals in SNP-PCA space, with locality labels:
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=as.factor(unlist(gl1$other)), #cex=size, pch=type of symbol, col=color
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(as.factor(unlist(gl1$other))), 
       pch=20, 
       col=as.factor(unique(unlist(gl1$other)))) 

# Perhaps we want to show disease status instead of locality:
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=as.factor(gl1$other$Trajectory), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(gl1$other$Trajectory), 
       pch=20, 
       col=as.factor(unique(gl1$other$Trajectory)))

# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(pca1$loadings[,1]),
            threshold=quantile(abs(pca1$loadings), 0.999)) 
##0.999 is the treshold we are observing, it can be another
threshold=quantile(abs(pca1$loadings), 0.999)  #0.1% threshold
threshold

# Get their locus names
gl1$loc.names[which(abs(pca1$loadings)>threshold)] 
# Only 6 surpassed the 0.1% threshold, NA (all under the line (see loading plot))

###DAPC: discriminant analysis on PCA, give a priori groups and ask how well the SNPS differentiate the groups
# Run the DAPC using disease status to group samples
disease.dapc <- dapc(gl1, pop=as.factor(unlist(gl1$other)), n.pca=8, n.da=3,
                     var.loadings=T, pca.info=T)
disease.dapc
# Scatterplot of results
scatter.dapc(disease.dapc, grp=as.factor(unlist(gl1$other)), legend=T)

# Plot the posterior assignment probabilities to each group
compoplot(disease.dapc)
disease.dapc

# Which loci contribute the most to distinguishing Healthy vs. Sick individuals?
loadingplot(abs(disease.dapc$var.load), 
            lab.jitter=1, 
            threshold=quantile(abs(disease.dapc$var.load), probs=0.999))
