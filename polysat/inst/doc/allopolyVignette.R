### R code from vignette source 'allopolyVignette.Rnw'

###################################################
### code chunk number 1: allopolyVignette.Rnw:125-130
###################################################
library(polysat)
data(AllopolyTutorialData)
summary(AllopolyTutorialData)
# make a copy of the dataset to modify
mydata <- AllopolyTutorialData


###################################################
### code chunk number 2: allopolyVignette.Rnw:147-156
###################################################
# Calculate the length of each genotype vector (= the number of alleles) and
# construct a TRUE/FALSE matrix of whether that number is greater than four.
tooManyAlleles <- apply(Genotypes(mydata), c(1,2), function(x) length(x[[1]])) > 4
# Find position(s) in the matrix that are TRUE.
which(tooManyAlleles, arr.ind=TRUE) # 43rd sample, second locus
# Look at the identified genotype, then replace it with missing data.
Genotype(mydata, 43, 2)
Genotype(mydata, 43, 2) <- Missing(mydata)
Genotype(mydata, 43, 2)


###################################################
### code chunk number 3: allopolyVignette.Rnw:163-165 (eval = FALSE)
###################################################
## mydist <- meandistance.matrix(mydata, distmetric=Lynch.distance,
##                               progress=FALSE)


###################################################
### code chunk number 4: allopolyVignette.Rnw:168-169
###################################################
load("vignettebuild/AllopolyTutorialDist.RData")


###################################################
### code chunk number 5: allopolyVignette.Rnw:172-175
###################################################
require(ape)
mynj <- nj(mydist)
plot(mynj, type="unrooted")


###################################################
### code chunk number 6: allopolyVignette.Rnw:185-186
###################################################
mydata <- deleteSamples(mydata, c("301","302","303"))


###################################################
### code chunk number 7: allopolyVignette.Rnw:192-196
###################################################
par(mfrow=c(2,1))
mypca <- cmdscale(mydist[Samples(mydata), Samples(mydata)])
plot(mypca[,1], mypca[,2])
hist(mypca[,1], breaks=30)


###################################################
### code chunk number 8: allopolyVignette.Rnw:204-206
###################################################
pop1ind <- Samples(mydata)[mypca[,1] <= 0]
pop2ind <- Samples(mydata)[mypca[,1] > 0]


###################################################
### code chunk number 9: allopolyVignette.Rnw:248-279
###################################################
nloc <- length(Loci(mydata)) # 7 loci
# lists to contain results of alleleCorrelations
CorrPop1 <- CorrPop2 <- list()
length(CorrPop1) <- length(CorrPop2) <- nloc
names(CorrPop1) <- names(CorrPop2) <- Loci(mydata)
# lists to contain results of testAlGroups
TAGpop1param1 <- list()
length(TAGpop1param1) <- nloc
names(TAGpop1param1) <- Loci(mydata)
TAGpop1param2 <- TAGpop1param3 <- TAGpop1param1
TAGpop2param1 <- TAGpop2param2  <- TAGpop2param3 <- TAGpop1param1

# loop through loci
for(L in Loci(mydata)){
  # allele correlations
  CorrPop1[[L]] <- alleleCorrelations(mydata, samples=pop1ind, locus = L)
  CorrPop2[[L]] <- alleleCorrelations(mydata, samples=pop2ind, locus = L)
  # default parameter set
  TAGpop1param1[[L]] <- testAlGroups(mydata, CorrPop1[[L]], samples=pop1ind)
  TAGpop2param1[[L]] <- testAlGroups(mydata, CorrPop2[[L]], samples=pop2ind)
  # optimized for homoplasy
  TAGpop1param2[[L]] <- testAlGroups(mydata, CorrPop1[[L]], samples=pop1ind,
                                     rare.al.check=0)
  TAGpop2param2[[L]] <- testAlGroups(mydata, CorrPop2[[L]], samples=pop2ind,
                                     rare.al.check=0)
  # optimized for null alleles
  TAGpop1param3[[L]] <- testAlGroups(mydata, CorrPop1[[L]], samples=pop1ind,
                                     null.weight=0)
  TAGpop2param3[[L]] <- testAlGroups(mydata, CorrPop2[[L]], samples=pop2ind,
                                     null.weight=0)
}


###################################################
### code chunk number 10: allopolyVignette.Rnw:291-293
###################################################
CorrPop1[["Loc6"]]$significant.pos
CorrPop2[["Loc6"]]$significant.pos


###################################################
### code chunk number 11: allopolyVignette.Rnw:303-304
###################################################
mydata <- deleteLoci(mydata, loci="Loc6")


###################################################
### code chunk number 12: allopolyVignette.Rnw:313-315
###################################################
# Population 1
heatmap(CorrPop1[["Loc1"]]$heatmap.dist, symm=TRUE)


###################################################
### code chunk number 13: allopolyVignette.Rnw:317-324
###################################################
# A plot to show how the colors correspond to p-values in the
# heat map; you can repeat this for the other heat maps in this
# tutorial if you wish.
plot(x=seq(min(CorrPop1[["Loc1"]]$heatmap.dist),
           max(CorrPop1[["Loc1"]]$heatmap.dist), length.out=12),
     y=rep(1,12), xlab="P-values", ylab="", bg=heat.colors(12), 
     pch=22, cex=3)


###################################################
### code chunk number 14: allopolyVignette.Rnw:326-331
###################################################
CorrPop1[["Loc1"]]$Kmeans.groups
CorrPop1[["Loc1"]]$UPGMA.groups
TAGpop1param1[["Loc1"]]$assignments
TAGpop1param2[["Loc1"]]$assignments
TAGpop1param3[["Loc1"]]$assignments


###################################################
### code chunk number 15: allopolyVignette.Rnw:333-335
###################################################
# Population 2
heatmap(CorrPop2[["Loc1"]]$heatmap.dist, symm=TRUE)


###################################################
### code chunk number 16: allopolyVignette.Rnw:337-342
###################################################
CorrPop2[["Loc1"]]$Kmeans.groups
CorrPop2[["Loc1"]]$UPGMA.groups
TAGpop2param1[["Loc1"]]$assignments
TAGpop2param2[["Loc1"]]$assignments
TAGpop2param3[["Loc1"]]$assignments


###################################################
### code chunk number 17: allopolyVignette.Rnw:355-357
###################################################
AssignToUse <- list()
AssignToUse[[1]] <- TAGpop1param2[["Loc1"]]


###################################################
### code chunk number 18: allopolyVignette.Rnw:364-419 (eval = FALSE)
###################################################
## heatmap(CorrPop1[["Loc2"]]$heatmap.dist, symm=TRUE)
## CorrPop1[["Loc2"]]$Kmeans.groups
## CorrPop1[["Loc2"]]$UPGMA.groups
## TAGpop1param1[["Loc2"]]$assignments
## TAGpop1param2[["Loc2"]]$assignments
## TAGpop1param3[["Loc2"]]$assignments
## # Population 2
## heatmap(CorrPop2[["Loc2"]]$heatmap.dist, symm=TRUE)
## CorrPop2[["Loc2"]]$Kmeans.groups
## CorrPop2[["Loc2"]]$UPGMA.groups
## TAGpop2param1[["Loc2"]]$assignments
## TAGpop2param2[["Loc2"]]$assignments
## TAGpop2param3[["Loc2"]]$assignments
## 
## heatmap(CorrPop1[["Loc3"]]$heatmap.dist, symm=TRUE)
## CorrPop1[["Loc3"]]$Kmeans.groups
## CorrPop1[["Loc3"]]$UPGMA.groups
## TAGpop1param1[["Loc3"]]$assignments
## TAGpop1param2[["Loc3"]]$assignments
## TAGpop1param3[["Loc3"]]$assignments
## # Population 2
## heatmap(CorrPop2[["Loc3"]]$heatmap.dist, symm=TRUE)
## CorrPop2[["Loc3"]]$Kmeans.groups
## CorrPop2[["Loc3"]]$UPGMA.groups
## TAGpop2param1[["Loc3"]]$assignments
## TAGpop2param2[["Loc3"]]$assignments
## TAGpop2param3[["Loc3"]]$assignments
## 
## heatmap(CorrPop1[["Loc4"]]$heatmap.dist, symm=TRUE)
## CorrPop1[["Loc4"]]$Kmeans.groups
## CorrPop1[["Loc4"]]$UPGMA.groups
## TAGpop1param1[["Loc4"]]$assignments
## TAGpop1param2[["Loc4"]]$assignments
## TAGpop1param3[["Loc4"]]$assignments
## # Population 2
## heatmap(CorrPop2[["Loc4"]]$heatmap.dist, symm=TRUE)
## CorrPop2[["Loc4"]]$Kmeans.groups
## CorrPop2[["Loc4"]]$UPGMA.groups
## TAGpop2param1[["Loc4"]]$assignments
## TAGpop2param2[["Loc4"]]$assignments
## TAGpop2param3[["Loc4"]]$assignments
## 
## heatmap(CorrPop1[["Loc5"]]$heatmap.dist, symm=TRUE)
## CorrPop1[["Loc5"]]$Kmeans.groups
## CorrPop1[["Loc5"]]$UPGMA.groups
## TAGpop1param1[["Loc5"]]$assignments
## TAGpop1param2[["Loc5"]]$assignments
## TAGpop1param3[["Loc5"]]$assignments
## # Population 2
## heatmap(CorrPop2[["Loc5"]]$heatmap.dist, symm=TRUE)
## CorrPop2[["Loc5"]]$Kmeans.groups
## CorrPop2[["Loc5"]]$UPGMA.groups
## TAGpop2param1[["Loc5"]]$assignments
## TAGpop2param2[["Loc5"]]$assignments
## TAGpop2param3[["Loc5"]]$assignments


###################################################
### code chunk number 19: allopolyVignette.Rnw:423-427
###################################################
AssignToUse[[2]] <- TAGpop1param1[["Loc2"]]
AssignToUse[[3]] <- TAGpop1param1[["Loc3"]]
AssignToUse[[4]] <- TAGpop1param1[["Loc4"]]
AssignToUse[[5]] <- TAGpop1param1[["Loc5"]]


###################################################
### code chunk number 20: allopolyVignette.Rnw:432-433
###################################################
heatmap(CorrPop1[["Loc7"]]$heatmap.dist, symm=TRUE)


###################################################
### code chunk number 21: allopolyVignette.Rnw:435-436
###################################################
heatmap(CorrPop2[["Loc7"]]$heatmap.dist, symm=TRUE)


###################################################
### code chunk number 22: allopolyVignette.Rnw:438-444
###################################################
TAGpop1param1[["Loc7"]]$assignments
TAGpop1param2[["Loc7"]]$assignments
TAGpop1param3[["Loc7"]]$assignments
TAGpop2param1[["Loc7"]]$assignments
TAGpop2param2[["Loc7"]]$assignments
TAGpop2param3[["Loc7"]]$assignments


###################################################
### code chunk number 23: allopolyVignette.Rnw:459-460
###################################################
mydata <- deleteLoci(mydata, loci="Loc7")


###################################################
### code chunk number 24: allopolyVignette.Rnw:472-475
###################################################
Loc1Param1Merged <- mergeAlleleAssignments(list(TAGpop1param1[["Loc1"]], 
                                                TAGpop2param1[["Loc1"]]))
Loc1Param1Merged[[1]]$assignments


###################################################
### code chunk number 25: allopolyVignette.Rnw:502-504
###################################################
recodedData <- recodeAllopoly(mydata, AssignToUse)
summary(recodedData)


###################################################
### code chunk number 26: allopolyVignette.Rnw:512-516
###################################################
for(L in Loci(recodedData)){
  proportionmissing <- mean(isMissing(recodedData, loci=L))
  cat(paste(L,":",proportionmissing,"missing"),sep="\n")
}


###################################################
### code chunk number 27: allopolyVignette.Rnw:524-525
###################################################
table(Ploidies(recodedData))


###################################################
### code chunk number 28: allopolyVignette.Rnw:541-549
###################################################
catResults <- list()
length(catResults) <- length(Loci(mydata))
names(catResults) <- Loci(mydata)

for(L in Loci(mydata)){
  cat(L, sep="\n")
  catResults[[L]] <- catalanAlleles(mydata, locus=L, verbose=TRUE)
}


