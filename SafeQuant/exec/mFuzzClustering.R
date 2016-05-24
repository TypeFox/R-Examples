# mFuzz cluster SafeQunat results
# 
# Author: erik ahrne
###############################################################################
### LOAD LIBRARIES
library(Mfuzz)
### LOAD LIBRARIES END

### LOAD DATA

### UPDATE THIS ###
sqRDataFile <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/UJenal/PabloManfredi/AdditionalFiles/Clustering/SQ_Results/SQ_Results_SQ.rData"
### UPDATE THIS ###
load(sqRDataFile)

eset <- sqaProtein$eset
eset <- eset[!fData(eset)$isFiltered,]

# remove all objects we don't need
rm(list=ls()[ls() != "eset"])

### LOAD DATA END

### PARAMS

### UPDATE THIS ###
pdfFile <- "//Volumes/pcf01$/Schmidt_Group/ProjectSQ/UJenal/PabloManfredi/AdditionalFiles/Clustering/mfuzz_clusters/clusters.pdf"
xlsExportFolder <- "/Volumes/pcf01$/Schmidt_Group/ProjectSQ/UJenal/PabloManfredi/AdditionalFiles/Clustering/mfuzz_clusters/"
### UPDATE THIS ###

nbClusters <- 8
min.mem <- .3

# graphics
colo <- "fancy" 
ylab <- "Standardized Intensity"
ylim <- c(-1.5,1.5)
xlab <- "Condition" 

### PARAMS END

# RE-LABEL, RE-ORDER ETC..

### expDesign overview
pData(eset)
#               condition isControl globalNormFactors
# 008_UJ2_SP3 Condition 1      TRUE         1.0000000
# 009_UJ2_SP3 Condition 2     FALSE         1.1006607
# 010_UJ2_SP3 Condition 3     FALSE         1.7372176
# 011_UJ2_SP3 Condition 4     FALSE         0.6141721
# 012_UJ2_SP3 Condition 5     FALSE         0.4752899
# 013_UJ2_SP3 Condition 6     FALSE         2.7911497

# MS1 intensity data overview
head(exprs(eset))
#                          008_UJ2_SP3 009_UJ2_SP3 010_UJ2_SP3 011_UJ2_SP3
# PA0001|ID:15869052|dnaA|   986552.56  1013825.06 1386022.385   807265.30
# PA0002|ID:15869053|dnaN|   600255.24   413466.61  586068.001   845047.10
# PA0004|ID:15869055|gyrB|  1461578.82  1162965.44  763840.304  1459000.42
# PA0006|ID:15871841|         21491.55    20931.19    7314.504    14283.51
# PA0007|ID:15869056|         30134.30    32335.09    5421.264    42139.59
# PA0008|ID:15871842|glyS|   627827.95   472942.78  379193.512   602734.84
#                          012_UJ2_SP3  013_UJ2_SP3
# PA0001|ID:15869052|dnaA|   654265.26 3711074.1043
# PA0002|ID:15869053|dnaN|   968620.05  458055.4280
# PA0004|ID:15869055|gyrB|  1473959.42 1162453.0830
# PA0006|ID:15871841|         14671.21    9532.0914
# PA0007|ID:15869056|         55099.15     619.1674
# PA0008|ID:15871842|glyS|   602605.97  378398.1889

# additional info overview 
# Note that the data is rolled up on the Protein level. I.e. the peptide column is a bit confusing.
# A protein is typically represented by multiple peptides (the listed peptide is the top-scoring peptide per protein)
str(fData(eset))
# 'data.frame':	2808 obs. of  17 variables:
#  $ proteinName       : Factor w/ 3161 levels "PA0001|ID:15869052|dnaA|",..: 1 2 3 5 6 7 8 9 10 11 ...
#  $ proteinDescription: Factor w/ 1550 levels "(2Fe-2S) ferredoxin [Pseudomonas aeruginosa PAO1]",..: 637 210 204 189 919 878 877 919 919 1173 ...
#  $ peptide           : Factor w/ 19436 levels "AAAAAFAPQLLDYK",..: 140 4192 16148 3009 2062 14587 10557 17359 7774 2663 ...
#  $ idScore           : num  690.7 1169.9 1917 41.6 81.8 ...
#  $ mass              : num  3410 1734 2248 1466 2835 ...
#  $ pMassError        : num  -0.621 -0.515 0.075 0.329 -0.023 ...
#  $ mz                : num  1138 868 1125 734 946 ...
#  $ retentionTime     : num  72.2 89.6 86.4 61.7 81.6 ...
#  $ charge            : num  3 2 2 2 3 2 2 2 2 2 ...
#  $ ptm               : Factor w/ 97 levels "","[10] Oxidation (M)",..: 1 1 1 1 1 1 1 1 1 1 ...
#  $ isNormAnchor      : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ isFiltered        : logi [1:2808(1d)] FALSE FALSE FALSE FALSE FALSE FALSE ...
#   ..- attr(*, "dimnames")=List of 1
#   .. ..$ : chr  "PA0001|ID:15869052|dnaA|" "PA0002|ID:15869053|dnaN|" "PA0004|ID:15869055|gyrB|" "PA0006|ID:15871841|" ...
#  $ nbPtmsPerPeptide  : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ allAccessions     : Factor w/ 3218 levels "PA0001|ID:15869052|dnaA|",..: 1 2 3 5 6 7 8 9 10 11 ...
#  $ nbPeptides        : int [1:2808(1d)] 10 14 21 1 3 14 3 2 1 2 ...
#   ..- attr(*, "dimnames")=List of 1
#   .. ..$ : chr  "PA0001|ID:15869052|dnaA|" "PA0002|ID:15869053|dnaN|" "PA0004|ID:15869055|gyrB|" "PA0006|ID:15871841|" ...
#  $ idQValue          : num  0 0 0 0.002326 0.000469 ...
#  $ allAccession      : chr  "PA0001|ID:15869052|dnaA|" "PA0002|ID:15869053|dnaN|" "PA0004|ID:15869055|gyrB|" "PA0006|ID:15871841|" ...

### RE-NAME SAMPLES  e.g. 008_UJ2_SP3 -> 008
rownames(pData(eset)) <- gsub("_UJ2.*","",rownames(pData(eset)))
# you can ofcourse name your sample whatever e.g.
#rownames(pData(eset))  <- c("s1","s2","s3","A","B","C")

# make sure expDesign is in agreement with exprs data
colnames(exprs(eset)) <- rownames(pData(eset))

# re-order 
pData(eset) <- pData(eset)[c("013","012","011","008","009","010"),]
# and make sure expDesign is in agreement with exprs data
exprs(eset) <- exprs(eset)[,rownames(pData(eset))]

# check
colnames(exprs(eset))
# [1] "013" "012" "011" "008" "009" "010"

# RE-LABEL, RE-ORDER ETC.. END





# log
exprs(eset) <- log10(exprs(eset))

# standardise
# good idea when including all data?
# e.g. only cluster data regulated at least X in at least one cond
eset <- standardise(eset)

### estimate fuzzification parameter using Schwammle and Jensen method
#       The algorithm needs a fuzzification parameter m in the range [1,n] which
#       determines the degree of fuzziness in the clusters. When
#       m reaches the value of 1 the algorithm works like a crisp
#       partitioning algorithm and for larger values of m the
#       overlapping of clusters is tend to be more.
m <- mestimate(eset)

### CLUSTER
#       This function is the core function for soft clustering. It groups genes based on the Euclidean distance
#       and the c-means objective function which is a weighted square error function. Each gene is assigned
#       a membership value between 0 and 1 for each cluster. Hence, genes can be assigned to different
#       clusters in a gradual manner. This contrasts hard clustering where each gene can belongs to a single
#       cluster.
#       Algorithm in brief:
#       1. Choose the number k of clusters.
#       2. Randomly generate k clusters and determine the cluster centers, or directly
#       generate k random points as cluster centers.
#       3. Assign each point to the nearest cluster center.
#       4. Recompute the new cluster centers.
#       5. Repeat the two previous steps until some convergence criterion is met (usually
#       that the assignment has not changed).
fuzzification <- mfuzz(eset,c=nbClusters,m=m)
### GRAPHICS
pdf(pdfFile)

### PLOT CLUSTERS WITH TREND LINE
par(mfrow=c(3,3))
for(clusterNb in 1:nbClusters){
	mfuzz.plot2(eset
					,cl=fuzzification
					, mfrow=NA
					, min.mem=min.mem
					, x11 = FALSE
					, ylim=ylim
					, single=clusterNb
					, ylab=ylab
					, xlab=xlab,colo=colo
					, time.labels=colnames(exprs(eset))
					,las=2)
			
	lines(1:ncol(eset),fuzzification$centers[clusterNb,], lwd=2)
}
mfuzzColorBar(col=colo,main="Membership",cex.main=1)
dev.off()
cat("Figures exported to", pdfFile, "\n")
### GRAPHICS END

### XLS EXPORT
for(i in 1:nbClusters){
	
	xlsExportFile = paste(xlsExportFolder,"cluster",i,".xls" ,sep="")
	protein <- names(fuzzification$cluster[fuzzification$cluster == i])
	description <- fData(eset)[protein,]$proteinDescription # get description of each protein
	
	membershipValue <- fuzzification$membership[fuzzification$cluster == i,i]
	outDf <- data.frame(protein=protein,description=description,membershipValue=as.vector(unlist(membershipValue)),exprs(eset)[protein,])
	outDf <- subset(outDf,membershipValue > min.mem)
	write.table(file=xlsExportFile,outDf, row.names=F,sep="\t")
	cat("Created file", xlsExportFile, "\n")
}
### XLS EXPORT END
cat("DONE\n")