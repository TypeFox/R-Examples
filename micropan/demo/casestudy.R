# This script contains all the code found in the casestudy.pdf document
# The comments are scarce, read the pdf-document for explanations.
#
# It was never our intention that this script should be run (sourced) as a whole,
# it will not work unless sourced in a directory set up exactly as described in 
# the casestudy document. Instead, use the code snippets below as examples that 
# you can copy and paste into your own R-script, and edit to fit your need.
#
# Lars Snipen & Kristian Hovde Liland



#--- Page 3 ---
# Reading the genome.table
library(micropan)
options(stringsAsFactors=FALSE)
genome.table <- read.table("data/Mpneumoniae.txt", sep="\t",
                           header=TRUE)



#--- Page 4 ---
# Downloading completed genomes from NCBI
for(i in 1:4){
  out.file <- file.path("data/genomes", genome.table$File[i])
  entrezDownload(genome.table$Accession[i], out.file)
}

# Downloading draft genomes from NCBI
for(i in 5:7){
  contig.acc <- getAccessions(genome.table$MasterRecord[i])
  out.file <- file.path("data/genomes", genome.table$File[i])
  entrezDownload(contig.acc,out.file)
}



#--- Page 5 ---
# Calling genes by Prodigal
for( i in 1:dim(genome.table)[1] ){
  cat("Predicting genes in", genome.table$File[i], "...\n")
  in.file <- file.path("data/genomes", genome.table$File[i])
  out.file <- file.path("data/proteins", genome.table$File[i])
  prodigalPredict(in.file, out.file)
}



#--- Page 6 ---
# Prepping protein files
for( i in 1:dim(genome.table)[1] ){
  cat("Preparing", genome.table$File[i], "...\n")
  in.file <- file.path("data/proteins", genome.table$File[i])
  gid <- genome.table$GID.tag[i]
  out.file <- file.path("data/prepped", genome.table$File[i])
  panPrep(in.file, gid, out.file)
}

# Reading FASTA file, displaying first 3 sequence headers
fdta <- readFasta("data/prepped/Mpneumoniae_M129_GID1.fsa")
fdta$Header[1:3]



#--- Page 7 ---
# Running BLAST
in.files <- file.path("data/prepped", dir("data/prepped"))
out.folder <- "blast"
blastAllAll(in.files, out.folder)



#--- Page 8 ---
# Running HMMER3
in.files <- file.path("data/prepped", dir("data/prepped"))
db <- "/usr/share/pfam/Pfam-A.hmm" # edit this to match your system
out.folder <- "pfam"
hmmerScan(in.files, db, out.folder)



#--- Page 9 ---
# Computing BLAST distances
blast.files <- file.path("blast", dir("blast"))
blast.distances <- bDist(blast.files)
save(blast.distances, file="res/blast_distances.RData")



#--- Page 10 ---
# Histogram of BLAST distances
hist(blast.distances[,3], breaks=50, col="tan4",
     xlab="BLAST distance", ylab="Number of distances")



#--- Page 11 ---
# Clustering based on BLAST distances
cluster.blast <- bClust(blast.distances, linkage="complete",
                        threshold=0.75)
# How many clusters?
length(unique(cluster.blast))



# --- Page 12 ---
# Display clustering of first 7 sequences
cluster.blast[1:7]



#--- Page 13 ---
# Reading results from hmmScan
pfam.files <- file.path("pfam", dir("pfam"))
pfam.table <- NULL
for(i in 1:length(pfam.files)){
  tab <- readHmmer(pfam.files[i])
  tab <- hmmerCleanOverlap(tab)
  pfam.table <- rbind(pfam.table, tab)
}
save(pfam.table, file="res/pfam_table.RData")

# Clustering based on domain sequences
cluster.pfam <- dClust(pfam.table)

# How many clusters?
length(unique(cluster.pfam))

# Display clustering of first 7 sequences
cluster.pfam[1:7]



#--- Page 14 ---
# cluster.info attribute of cluster 1
attr(cluster.pfam, "cluster.info")[1]



#--- Page 15 ---
# Constructing pan-matrices
pm.blast <- panMatrix(cluster.blast)
pm.pfam <- panMatrix(cluster.pfam)

# Plotting pan-matrices
par(mfrow=c(2,1))
plot(pm.blast)
plot(pm.pfam)

# Summary of pan-matrices
summary(pm.blast)
summary(pm.pfam)



#--- Page 16 ---
# Creating pangenome trees
blast.tree <- panTree(pm.blast, dist.FUN=distManhattan)
pfam.tree <- panTree(pm.pfam)

# Plotting the trees
par(mfrow=c(1,2))
plot(blast.tree)
plot(pfam.tree)



#--- Page 18 ---
# Nicer pangenome tree
blast.tree <- panTree(pm.blast, nboot=100) # tree with bootstrapping
my.lab <- genome.table$Strain
names( my.lab ) <- genome.table$GID.tag
my.col <- genome.table$Color
names( my.col ) <- genome.table$GID.tag
plot(blast.tree, leaf.lab=my.lab, col=my.col,
     xlab="Manhattan distances",
     main="Pangenome tree for Mycoplasma pneumoniae")



#--- Page 19 ---
# Fitting binomial mixture model
binomix <- binomixEstimate(pm.blast, K.range=2:7)

# Displaying the BIC.table
binomix$BIC.table

# Summary of model
summary(binomix)

# Chao estimate
pan.size <- chao(pm.blast)
pan.size



#--- Page 20 ---
# Plotting binomial mixture
plot(binomix)



#--- Page 21 ---
# Computing rarefaction and plotting
r <- rarefaction(pm.blast, n.perm=100)
plot(r)



#--- Page 22 ---
# Fitting Heaps model
h <- heaps(pm.blast, n.perm=100)
h

# Computing fluidity
f <- fluidity(pm.blast, n.sim=100)
f

# Computing Jaccard distances
J <- distJaccard(pm.blast)
mean(J)



#--- Page 23 ---
# Histogram of Jaccard distances
hist(J, breaks=10, col="tan")




