###############################################################################
## Illumina's default method vs. robloxbioc for technical replicates
###############################################################################

###############################################################################
## References:
## Dunning, M.J., Smith, M.L., Ritchie, M.E., Tavare, S.:
## beadarray: R classes and methods for Illumina bead-based data. 
## Bioinformatics 2007, 23(16):2183-4.
##
## M.J. Dunning, N.L. Barbosa-Morais, A.G. Lynch, S. Tavaré and M.E. Ritchie:
## Statistical issues in the analysis of Illumina data.
## BMC Bioinformatics 2008, 9:85.
###############################################################################

###############################################################################
## Data:
## Can be obtained via
## http://www.compbio.group.cam.ac.uk/Resources/spike/index.html
###############################################################################

## Load the required packages
library(beadarray)
library(RobLoxBioC)

###############################################################################
## Extract all *.zip file to directory "SpikeInData".
## Copy spike_targets.txt to directory "SpikeInData".
##
## Code to read the bead level data from the directory "SpikeInData"
##
## NOTE: reading in the raw data for the entire experiment requires at
## least 4Gb of RAM for each processing method.  
###############################################################################

###########################################################
## Read targets information
targets <- read.table("./SpikeInData/spike_targets.txt",header=TRUE)
arraynms <- as.character(targets$ArrayNo)

## Use sharpened, subtracted data from text files
spikeInData <- readIllumina(path = "./SpikeInData", arrayNames=arraynms[1:2], 
                            useImages=FALSE, textType=".csv")
#save(spikeInData, compress = TRUE, file = "spikeInData.RData")
#load(file = "spikeInData.RData")

## takes about 80 sec (Core i5 520M with 8 GByte RAM)
system.time(res.ill <- createBeadSummaryData(spikeInData, log = TRUE, imagesPerArray = 2))
#save(res.ill, file = "SpikeInData_Ill.RData")
#load(file = "SpikeInData_Ill.RData")

## takes about 490 sec (Core i5 520M with 8 GByte RAM)
system.time(res.rmx <- robloxbioc(spikeInData, imagesPerArray = 2))
#save(res.rmx, file = "SpikeInData_rmx.RData")
#load(file = "SpikeInData_rmx.RData")

###########################################################
## From Dunning et al. (2008):
## "The spikes were added at concentrations of 1000, 300, 100, 30, 10 and 3 pM 
## on the six arrays from the first four BeadChips. A further four chips were 
## hybridised with spikes at concentrations of 1, 0.3, 0.1, 0.03, 0.01 and 0 pM. 
## The spikes on a given array were all added at the same concentration. Each 
## concentration was allocated to the same position on all replicate BeadChips. 
## For example, 1000 pM was always array 1 on a chip and 300 pM was array 2 and 
## so on."
###########################################################

A <- c(1,7,13,19)
cor.ill.A <- cor(exprs(res.ill[,A]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.A <- cor(exprs(res.rmx[,A]), method = "spearman", use = "pairwise.complete.obs")
(diff.A <- cor.rmx.A-cor.ill.A)
(rel.A <- cor.rmx.A/cor.ill.A)
range(diff.A[col(diff.A) > row(diff.A)])
100*(range(rel.A[col(rel.A) > row(rel.A)])-1)

cor.ill.A1 <- cor(exprs(res.ill[,A]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.A1 <- cor(exprs(res.rmx[,A]), method = "pearson", use = "pairwise.complete.obs")
(diff.A1 <- cor.rmx.A1-cor.ill.A1)
(rel.A1 <- cor.rmx.A1/cor.ill.A1)
range(diff.A1[col(diff.A1) > row(diff.A1)])
100*(range(rel.A1[col(rel.A1) > row(rel.A1)])-1)


B <- c(2,8,14,20)
cor.ill.B <- cor(exprs(res.ill[,B]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.B <- cor(exprs(res.rmx[,B]), method = "spearman", use = "pairwise.complete.obs")
(diff.B <- cor.rmx.B-cor.ill.B)
(rel.B <- cor.rmx.B/cor.ill.B)
range(diff.B[col(diff.B) > row(diff.B)])
100*(range(rel.B[col(rel.B) > row(rel.B)])-1)

cor.ill.B1 <- cor(exprs(res.ill[,B]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.B1 <- cor(exprs(res.rmx[,B]), method = "pearson", use = "pairwise.complete.obs")
(diff.B1 <- cor.rmx.B1-cor.ill.B1)
(rel.B1 <- cor.rmx.B1/cor.ill.B1)
range(diff.B1[col(diff.B1) > row(diff.B1)])
100*(range(rel.B1[col(rel.B1) > row(rel.B1)])-1)


C <- c(3,9,15,21)
cor.ill.C <- cor(exprs(res.ill[,C]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.C <- cor(exprs(res.rmx[,C]), method = "spearman", use = "pairwise.complete.obs")
(diff.C <- cor.rmx.C-cor.ill.C)
(rel.C <- cor.rmx.C/cor.ill.C)
range(diff.C[col(diff.C) > row(diff.C)])
100*(range(rel.C[col(rel.C) > row(rel.C)])-1)

cor.ill.C1 <- cor(exprs(res.ill[,C]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.C1 <- cor(exprs(res.rmx[,C]), method = "pearson", use = "pairwise.complete.obs")
(diff.C1 <- cor.rmx.C1-cor.ill.C1)
(rel.C1 <- cor.rmx.C1/cor.ill.C1)
range(diff.C1[col(diff.C1) > row(diff.C1)])
100*(range(rel.C1[col(rel.C1) > row(rel.C1)])-1)


D <- c(4,10,16,22)
cor.ill.D <- cor(exprs(res.ill[,D]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.D <- cor(exprs(res.rmx[,D]), method = "spearman", use = "pairwise.complete.obs")
(diff.D <- cor.rmx.D-cor.ill.D)
(rel.D <- cor.rmx.D/cor.ill.D)
range(diff.D[col(diff.D) > row(diff.D)])
100*(range(rel.D[col(rel.D) > row(rel.D)])-1)

cor.ill.D1 <- cor(exprs(res.ill[,D]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.D1 <- cor(exprs(res.rmx[,D]), method = "pearson", use = "pairwise.complete.obs")
(diff.D1 <- cor.rmx.D1-cor.ill.D1)
(rel.D1 <- cor.rmx.D1/cor.ill.D1)
range(diff.D1[col(diff.D1) > row(diff.D1)])
100*(range(rel.D1[col(rel.D1) > row(rel.D1)])-1)


E <- c(5,11,17,23)
cor.ill.E <- cor(exprs(res.ill[,E]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.E <- cor(exprs(res.rmx[,E]), method = "spearman", use = "pairwise.complete.obs")
(diff.E <- cor.rmx.E-cor.ill.E)
(rel.E <- cor.rmx.E/cor.ill.E)
range(diff.E[col(diff.E) > row(diff.E)])
100*(range(rel.E[col(rel.E) > row(rel.E)])-1)

cor.ill.E1 <- cor(exprs(res.ill[,E]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.E1 <- cor(exprs(res.rmx[,E]), method = "pearson", use = "pairwise.complete.obs")
(diff.E1 <- cor.rmx.E1-cor.ill.E1)
(rel.E1 <- cor.rmx.E1/cor.ill.E1)
range(diff.E1[col(diff.E1) > row(diff.E1)])
100*(range(rel.E1[col(rel.E1) > row(rel.E1)])-1)


F <- c(6,12,18,24)
cor.ill.F <- cor(exprs(res.ill[,F]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.F <- cor(exprs(res.rmx[,F]), method = "spearman", use = "pairwise.complete.obs")
(diff.F <- cor.rmx.F-cor.ill.F)
(rel.F <- cor.rmx.F/cor.ill.F)
range(diff.F[col(diff.F) > row(diff.F)])
100*(range(rel.F[col(rel.F) > row(rel.F)])-1)

cor.ill.F1 <- cor(exprs(res.ill[,F]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.F1 <- cor(exprs(res.rmx[,F]), method = "pearson", use = "pairwise.complete.obs")
(diff.F1 <- cor.rmx.F1-cor.ill.F1)
(rel.F1 <- cor.rmx.F1/cor.ill.F1)
range(diff.F1[col(diff.F1) > row(diff.F1)])
100*(range(rel.F1[col(rel.F1) > row(rel.F1)])-1)


#######################################
## second series of dilutions
#######################################
A <- c(25,31,37,43)
cor.ill.A <- cor(exprs(res.ill[,A]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.A <- cor(exprs(res.rmx[,A]), method = "spearman", use = "pairwise.complete.obs")
(diff.A <- cor.rmx.A-cor.ill.A)
(rel.A <- cor.rmx.A/cor.ill.A)
range(diff.A[col(diff.A) > row(diff.A)])
100*(range(rel.A[col(rel.A) > row(rel.A)])-1)

cor.ill.A1 <- cor(exprs(res.ill[,A]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.A1 <- cor(exprs(res.rmx[,A]), method = "pearson", use = "pairwise.complete.obs")
(diff.A1 <- cor.rmx.A1-cor.ill.A1)
(rel.A1 <- cor.rmx.A1/cor.ill.A1)
range(diff.A1[col(diff.A1) > row(diff.A1)])
100*(range(rel.A1[col(rel.A1) > row(rel.A1)])-1)


B <- c(26,32,38,44)
cor.ill.B <- cor(exprs(res.ill[,B]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.B <- cor(exprs(res.rmx[,B]), method = "spearman", use = "pairwise.complete.obs")
(diff.B <- cor.rmx.B-cor.ill.B)
(rel.B <- cor.rmx.B/cor.ill.B)
range(diff.B[col(diff.B) > row(diff.B)])
100*(range(rel.B[col(rel.B) > row(rel.B)])-1)

cor.ill.B1 <- cor(exprs(res.ill[,B]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.B1 <- cor(exprs(res.rmx[,B]), method = "pearson", use = "pairwise.complete.obs")
(diff.B1 <- cor.rmx.B1-cor.ill.B1)
(rel.B1 <- cor.rmx.B1/cor.ill.B1)
range(diff.B1[col(diff.B1) > row(diff.B1)])
100*(range(rel.B1[col(rel.B1) > row(rel.B1)])-1)


C <- c(27,33,39,45)
cor.ill.C <- cor(exprs(res.ill[,C]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.C <- cor(exprs(res.rmx[,C]), method = "spearman", use = "pairwise.complete.obs")
(diff.C <- cor.rmx.C-cor.ill.C)
(rel.C <- cor.rmx.C/cor.ill.C)
range(diff.C[col(diff.C) > row(diff.C)])
100*(range(rel.C[col(rel.C) > row(rel.C)])-1)
100*(1/rel.C[col(rel.C) > row(rel.C)]-1)


cor.ill.C1 <- cor(exprs(res.ill[,C]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.C1 <- cor(exprs(res.rmx[,C]), method = "pearson", use = "pairwise.complete.obs")
(diff.C1 <- cor.rmx.C1-cor.ill.C1)
(rel.C1 <- cor.rmx.C1/cor.ill.C1)
range(diff.C1[col(diff.C1) > row(diff.C1)])
100*(range(rel.C1[col(rel.C1) > row(rel.C1)])-1)
100*(1/rel.C1[col(rel.C1) > row(rel.C1)]-1)


D <- c(28,34,40,46)
cor.ill.D <- cor(exprs(res.ill[,D]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.D <- cor(exprs(res.rmx[,D]), method = "spearman", use = "pairwise.complete.obs")
(diff.D <- cor.rmx.D-cor.ill.D)
(rel.D <- cor.rmx.D/cor.ill.D)
range(diff.D[col(diff.D) > row(diff.D)])
100*(range(rel.D[col(rel.D) > row(rel.D)])-1)
100*(1/rel.D[col(rel.D) > row(rel.D)]-1)

cor.ill.D1 <- cor(exprs(res.ill[,D]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.D1 <- cor(exprs(res.rmx[,D]), method = "pearson", use = "pairwise.complete.obs")
(diff.D1 <- cor.rmx.D1-cor.ill.D1)
(rel.D1 <- cor.rmx.D1/cor.ill.D1)
range(diff.D1[col(diff.D1) > row(diff.D1)])
100*(range(rel.D1[col(rel.D1) > row(rel.D1)])-1)
100*(1/rel.D1[col(rel.D1) > row(rel.D1)]-1)


E <- c(29,35,41,47)
cor.ill.E <- cor(exprs(res.ill[,E]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.E <- cor(exprs(res.rmx[,E]), method = "spearman", use = "pairwise.complete.obs")
(diff.E <- cor.rmx.E-cor.ill.E)
(rel.E <- cor.rmx.E/cor.ill.E)
range(diff.E[col(diff.E) > row(diff.E)])
100*(range(rel.E[col(rel.E) > row(rel.E)])-1)

cor.ill.E1 <- cor(exprs(res.ill[,E]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.E1 <- cor(exprs(res.rmx[,E]), method = "pearson", use = "pairwise.complete.obs")
(diff.E1 <- cor.rmx.E1-cor.ill.E1)
(rel.E1 <- cor.rmx.E1/cor.ill.E1)
range(diff.E1[col(diff.E1) > row(diff.E1)])
100*(range(rel.E1[col(rel.E1) > row(rel.E1)])-1)


F <- c(30,36,42,48)
cor.ill.F <- cor(exprs(res.ill[,F]), method = "spearman", use = "pairwise.complete.obs")
cor.rmx.F <- cor(exprs(res.rmx[,F]), method = "spearman", use = "pairwise.complete.obs")
(diff.F <- cor.rmx.F-cor.ill.F)
(rel.F <- cor.rmx.F/cor.ill.F)
range(diff.F[col(diff.F) > row(diff.F)])
100*(range(rel.F[col(rel.F) > row(rel.F)])-1)

cor.ill.F1 <- cor(exprs(res.ill[,F]), method = "pearson", use = "pairwise.complete.obs")
cor.rmx.F1 <- cor(exprs(res.rmx[,F]), method = "pearson", use = "pairwise.complete.obs")
(diff.F1 <- cor.rmx.F1-cor.ill.F1)
(rel.F1 <- cor.rmx.F1/cor.ill.F1)
range(diff.F1[col(diff.F1) > row(diff.F1)])
100*(range(rel.F1[col(rel.F1) > row(rel.F1)])-1)



