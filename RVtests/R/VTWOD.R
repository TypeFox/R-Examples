VTWOD <-
function(x, y, polyphen.weight, flipPhenotype = 0, 
	npermutation = 1000, npermutation.max, min.nonsignificant.counts)
{   
#####################################################
# This script implements the VT test for pooled association of rare variants with a phenotype.
# See Price et al. AJHG 2010
# The script also includes implementations of T1, T5, and WE (Madsen-Browning) 
#     tests, optionally weighted with PolyPhen scores.
# 2011 - included new test designed for two-tailed variants (WOD at 1% and 5%)
#
# For speed, the script supports three modes of running: local on single CPU, 
#     multicore, and cluster (see options below).
#
# NOTE: currently, the script is configured to run the tests on a single gene.
# 
#   <npermutation> is an integer number of npermutation to perform
#   <phenotypeFile> is the name of a file with lines of format: individualID phenotypeValue
#   <snpWeightFile> is the name of a file with lines of format (weight is between 0 and 1): snpid weight
#   <genotypeFile> is the name of a file with lines of format (genotype is the number of rare alleles of the snp in the individual, typically one of {0,1,2}): individualID snpid genotype
#   <seed> is an optional random seed value
# 
# Output: p-values for
# score1, score1P - test T1 and T1P (see paper)
# score2, score2P - test T5 and T5P (see paper)
# score3, score3P - test WE and WEP (see paper)
# score4, score4P - test VT and VTP (see paper)
# WOD.01 - WOD method with 1% threshold
# WOD.05 - WOD method with 5% threshold
#
# This R implementation by Adam Kiezun, based on reference implementation in C by Alkes Price.
# Added WOD tests to the program in 2011 by Celia Greenwood
#######################################################
#
#source("E:/Changjiang Xu/raresnps/WOD_program/WOD.VariantTestsCJ.R")
#
#indFile <- "sample.pheno.txt"  # indid, pheno 
#snpFile <- "sample.polyphen2.weight"       # polyphen scores-set to 0.5 if missing 
#genoFile <- "sample.geno"         # file with genotypes
#                                         # indid, snpid, count(0,1,2) rare alleles
#flipPhenotype <- 0      # logical, should the phenotype be multiplied by -1
#npermutation <- 1000    # number of npermutation to perform
#seed <-  10985          # random seed, integer
#set.seed(seed)
#
# read data and remove problem individuals
#DIR<- "E:/Changjiang Xu/raresnps/WOD_program/"
#ind <- read.table(paste(DIR, indFile, sep=""), col.names=c("indid", "pheno"))
#csnp <- read.table(paste(DIR, snpFile, sep=""), col.names=c("snpid", "polyphen"))
#cgeno <- read.table(paste(DIR, genoFile, sep=""), col.names=c("indid", "snpid", "count"))
#source("E:/Changjiang Xu/raresnps/WOD_program/count2geno.R")
#source("E:/Changjiang Xu/raresnps/WOD_program/geno2count.R")
#geno<- count2geno(cgeno)
#cgeno<- geno2count(geno)
#
#######################################################
#Data formatting
if (is.null(colnames(x))) colnames(x)<- paste("snp", 1:ncol(x), sep="")
if (is.null(rownames(x))) rownames(x)<- 1:nrow(x)
cgeno <- geno2count(x)
ind <- data.frame(indid = rownames(x), pheno = y, stringsAsFactors = FALSE)
if (missing(polyphen.weight)) csnp<- data.frame(snpid = colnames(x), polyphen = rep(1, ncol(x)), stringsAsFactors = FALSE)
else {csnp<- as.data.frame(polyphen.weight, stringsAsFactors = FALSE)
	colnames(csnp)<- c("snpid", "polyphen")
	}

#permutation parameters
if (missing(min.nonsignificant.counts)) min.nonsignificant.counts<- 10  #?10, ?20
if (missing(npermutation.max)) npermutation.max<- npermutation


#######################################################
#The following codes are same with WOD.VariantTests.R except adding adaptive permutation.
#######################################################



#Sometimes phenotypes are annotated the opposite way of what we're expecting. If yes, then flip.
ind <- ind[,c("indid", "pheno")]
if (flipPhenotype){
  cat("flipping phenotypes\n")
  ind$pheno <- -1.0*ind$pheno
}
meanpheno <- mean(ind$pheno)

# VT code
#For each SNP, how many times it is seen
csnp$counts <- sapply(csnp$snpid, function(x){ sum(cgeno[cgeno$snpid == x,]$count) })

#For a SNP with total count c, how many counts are lower than c? (ie what is the rank or c in the order of counts)
csnp$countg <- sapply(csnp$counts, function(x){ length(unique(csnp[csnp$counts < x,]$counts)) })

#Sample size
N <- dim(ind)[1]

m1 <- merge(cgeno, csnp, by=c("snpid"))

#adjust polyphen scores
if (length(m1[(m1$counts >= N/50) & (m1$polyphen < 1.0),]$polyphen) > 0){
 m1[(m1$counts >= N/50) & (m1$polyphen < 1.0),]$polyphen <- 0.5
}

#pre-compute metrics that are independent of npermutation
m1$countSquare <- m1$count*m1$count
m1$countPolyphen <- m1$count*m1$polyphen
m1$countSquarePolyphenSquare <- m1$countPolyphen*m1$countPolyphen
f <- (1+m1$counts)/(2+2*N)
m1$weight <- 1/sqrt(f*(1.0-f))
m1$countWeight <- m1$count*m1$weight

#Create a single table by joining SNPs and genotypes by SNPid, and joining individuals by individual ID
m <- merge(m1, ind, by=c("indid"))
m <- m[m$counts < N,] #ignore common variants

#Compute sum for subsets of indices (the subsets are pre-computed)
mysum <- function(X, range, arr, whiches){
  for (i in range) { arr[i] <- sum(X[whiches[[i]]])}
  arr
}

#To improve speed, pre-compute everything that is independent of npermutation
ctg <- m$countg
mCount <- m$count;
mPolyphen <- m$polyphen;
mCountWeight <- m$countWeight
indPheno <- ind$pheno

nx <- length(ctg)
fctg <- as.factor(list(ctg)[[1]])
index <- fctg
one <- 1L
group <- rep.int(one, nx) + one * (as.integer(index) - one)
len <- length(unique(group))
arr <- double(len)
oneToLen <- 1:len
Msize <- dim(m)[1]

whiches <- vector("list", len)
for (i in oneToLen) { whiches[[i]] <- which(group == i) }

count <- mysum(mCount, range=oneToLen, arr=arr, whiches=whiches)
countSquare <- mysum(m$countSquare, range=oneToLen, arr=arr, whiches=whiches)
countSquarePolyphenSquare <- mysum(m$countSquarePolyphenSquare, range=oneToLen, arr=arr, whiches=whiches)
countPolyphen <- mysum(m$countPolyphen, range=oneToLen, arr=arr, whiches=whiches)

#Indices of variants below frequency thresholds
mBelow50 <- which(m$counts < N/50)
mBelow10 <- which(m$counts < N/10)

#Pre-compute cumulative sums, for VT test
csCount <- cumsum(count)
csCountSquare <- cumsum(countSquare)
csCountPolyphenMeanpheno <- cumsum(countPolyphen * meanpheno)
csCountSquarePolyphenSquare <- cumsum(countSquarePolyphenSquare)
csCountMeanpheno <- csCount*meanpheno
sqrtCsCountSquare <- sqrt(csCountSquare)
sqrtCsCountSquarePolyphenSquare <- sqrt(csCountSquarePolyphenSquare)

#Indices of individuals from m in ind (may be duplicate)
matchIds <- match(m$indid, ind$indid)

# Items needed for CG's implementation of Liu & Leal
   csnp$rare01 = (csnp$counts/(2*N)<0.01)  # 1% threshold for rare
   csnp$rare05 = (csnp$counts/(2*N)<0.05)  # 5% threshold for rare
   # find all individuals who carry no rare variants at 1% or 5%
   ind$norare01 = sapply(ind$indid,  function(x)  {
        all(is.na(match(m1[m1$indid==x,]$snpid, csnp$snpid[csnp$rare01] )))})
   ind$norare05 = sapply(ind$indid,  function(x)  {
        all(is.na(match(m1[m1$indid==x,]$snpid, csnp$snpid[csnp$rare05] )))})

################################################################
#Compute the test scores for many tests, for 1 permutation
getScores <- function(permute){
  if (permute){
    pheno.cg <- sample(ind$pheno)
    pheno <- pheno.cg[matchIds]
  } else {
    pheno.cg <- ind$pheno
    pheno <- pheno.cg[matchIds]
  }
  phenoCount         <- pheno * mCount
  phenoCountPolyphen <- phenoCount * mPolyphen
  phenoCountWeight   <- pheno * mCountWeight

  #Scores that count only rare variants, optionally weighted
  score1 <- sum(phenoCount[mBelow50])
  score2 <- sum(phenoCount[mBelow10])
  score1P <- sum(phenoCountPolyphen[mBelow50])
  score2P <- sum(phenoCountPolyphen[mBelow10])

  #Madsen-Browning score, optionally weighted
  score3 <- sum(phenoCountWeight)
  score3P <- sum(phenoCountWeight * mPolyphen)

  #VT test, optionally weighted
  #Aggregate for each count, to find the optimal threshold for VT test
  csPhenoCount <- cumsum(mysum(phenoCount, range=oneToLen, arr=arr, whiches=whiches))
  csPhenoCountPolyphen <- cumsum(mysum(phenoCountPolyphen, range=oneToLen, arr=arr, whiches=whiches))
  
  score4 <- max((csPhenoCount-csCountMeanpheno)/sqrtCsCountSquare)
  score4P <- max((csPhenoCountPolyphen-csCountPolyphenMeanpheno)/sqrtCsCountSquarePolyphenSquare)

  # new code by CG for Liu & Leal extension
  pheno.mean <- mean(pheno.cg[ind$norare01])
  pheno.sd <- sd(pheno.cg[ind$norare01])
  quant <- pnorm(pheno.cg, pheno.mean, pheno.sd)
  probnonnull <- 1 - 2*(quant*(quant<0.5) + (1-quant)*(quant>=0.5))
  stat.b <- t(abs(pheno.cg-meanpheno)) %*% probnonnull / sum(probnonnull)
  WOD.01 <- -stat.b  # do this because score expected to get larger
                     # under permuted null.  
  pheno.mean <- mean(pheno.cg[ind$norare05])
  pheno.sd <- sd(pheno.cg[ind$norare05])
  quant <- pnorm(pheno.cg, pheno.mean, pheno.sd)
  probnonnull <- 1 - 2*(quant*(quant<0.5) + (1-quant)*(quant>=0.5))
  stat.b <- t(abs(pheno.cg-meanpheno)) %*% probnonnull / sum(probnonnull)
  WOD.05 <- -stat.b    

  c(score1=score1, score1P=score1P, score2=score2, score2P=score2P,
    score3=score3, score3P=score3P, score4=score4, 
    score4P=score4P,  WOD.01 = WOD.01, WOD.05 = WOD.05)
}

##################################################################
#For a specific score, returns how often permuted data has a higher score than unpermuted data.
permwins <- function(scores.df, scorename) {
  unpermuted <- (unpermutedScores[c(scorename)])[[1]]
  ceiling(sum(scores.df[,c(scorename)] > unpermuted) + 0.5*sum(scores.df[,c(scorename)] == unpermuted))
}

##################################################################
#For all scores, returns how often permuted data has a higher score than unpermuted data.
PermLoop <- function(npermutation) {
	countwin<- 0
	pm<- 0
	while ((pm < npermutation) | ((pm < npermutation.max) & (min(countwin) < min.nonsignificant.counts)))
	{
	pm<- pm + 1
        scores <- getScores(permute=TRUE)
        scores.df <- as.data.frame(t(scores))
        pw1  <- permwins(scores.df, "score1")
        pw1P <- permwins(scores.df, "score1P")
        pw2  <- permwins(scores.df, "score2")
        pw2P <- permwins(scores.df, "score2P")
        pw3  <- permwins(scores.df, "score3")
        pw3P <- permwins(scores.df, "score3P")
        pw4  <- permwins(scores.df, "score4")
        pw4P <- permwins(scores.df, "score4P")
        pw.nll01 <- permwins(scores.df, "WOD.01")
        pw.nll05 <- permwins(scores.df, "WOD.05")
    if (pm==1) {
      permresult <- 
        data.frame(t( c(score1=pw1, score1P=pw1P, score2=pw2, 
          score2P=pw2P, score3=pw3, 
          score3P=pw3P, score4=pw4, score4P=pw4P, WOD.01=pw.nll01,
          WOD.05=pw.nll05))) }  else {
      permresult <- rbind(permresult,
          data.frame(t( c(score1=pw1, score1P=pw1P, score2=pw2, 
          score2P=pw2P, score3=pw3, 
          score3P=pw3P, score4=pw4, score4P=pw4P, WOD.01=pw.nll01,
          WOD.05=pw.nll05)))) }

	countwin<- colSums(permresult)
 }
 list(permresult = permresult, pm =pm)
}
############################################################

#P-values
pval <- function(permwins){ (permwins+1)/(npermutation+1) }

#################################################################

####################################
# Function Calls start here
##################################
#Unpermuted data for which we're looking for pvalues
unpermutedScores <- as.data.frame(t(getScores(permute=FALSE)))
#print(unpermutedScores)

# Run npermutation
pwpm<- PermLoop(npermutation)
pw <- pwpm$permresult
totalpm<- pwpm$pm
pw.df <- as.data.frame(t(apply(pw, 2, sum))) 

#cat("counts where permuted statistics are larger than unpermuted data\n")
#print(pw.df)
#cat("p-values\n")

#options(max.levels = NULL)
#options(width = 500)
#print(pval(pw.df), max.levels = NULL, width = 500)

list(score = unpermutedScores, nonsignificant.counts = pw.df, pvalue.empirical = pval(pw.df), pvalue.nominal = NA, total.permutation = totalpm)
}

