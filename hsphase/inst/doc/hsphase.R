### R code from vignette source 'hsphase.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: imageplot
###################################################

set.seed(718)
library(hsphase)
halfsibs <- .simulateHalfsib(numInd = 20)
imageplot(bmh(halfsibs),title = "Imageplot of simulated half-sib family")


###################################################
### code chunk number 2: rplot
###################################################
rplot(halfsibs,sort(sample(c(1:(1000*ncol(halfsibs))),size=ncol(halfsibs))))
title("Recombination of simulated half-sib family")


###################################################
### code chunk number 3: heatmap
###################################################
a <- .simulateHalfsib()
b <- .simulateHalfsib()
d <- rbind(a,b)
library(hsphase)
oh <- ohg(d)
heatmap(oh,symm=T,col=gray.colors(16,start=0,end=1),RowSideColors=as.character(c(rep(1,40),rep(2,40))),ColSideColors=as.character(c(rep(1,5),rep(2,40),rep(1,35))))


###################################################
### code chunk number 4: ohplot
###################################################
set.seed(100)
chr <- list()
sire <- list()
set.seed(1)
chr <- list()
for(i in 1:5)
{
	chr[[i]] <- .simulateHalfsib(numInd = 20, numSNP = 5000, recbound = 1:10)
	sire[[i]] <- ssp(bmh(chr[[i]]),chr[[i]])
	sire[[i]] <- sire[[i]][1,]+sire[[i]][2,]
	sire[[i]][sire[[i]]==18] <- 9
}

Genotype <- do.call(rbind, chr)
rownames(Genotype) <- 6:(nrow(Genotype)+5)
sire <- do.call(rbind, sire)
rownames(sire) <- 1:5
Genotype <- rbind(sire, Genotype)
oh <- ohg(Genotype)  # creating the Opposing Homozygote matrix
pedigree <- as.matrix(data.frame( c(1:5,6:(nrow(Genotype))),rep = c(rep(0,5), rep(1:5,rep(20,5)))))
ohplot(oh, Genotype, pedigree, check = TRUE)


