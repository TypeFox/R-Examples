### R code from vignette source 'sdcTable.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sdcTable.Rnw:87-88
###################################################
library(sdcTable)


###################################################
### code chunk number 2: sdcTable.Rnw:92-103
###################################################
# daten laden
workDat <- get(load('workDat.RData'))
inputData <- workDat$inputData
microData <- inputData$microDat
aggregatedData <- inputData$aggDat
dimList <- inputData$dimList

resOPT <- workDat$resOPT

completeData <- getInfo(workDat$resOPT, 'finalData')
completeData <- completeData[,-ncol(completeData)]


###################################################
### code chunk number 3: sdcTable.Rnw:124-125
###################################################
print(head(microData), row.names=FALSE)


###################################################
### code chunk number 4: sdcTable.Rnw:141-144
###################################################
lev.V1 <- as.character(sort(unique(microData$V1)))
lev.V2 <- as.character(sort(unique(microData$V2)))
lev.V3 <- as.character(sort(unique(microData$V3)))


###################################################
### code chunk number 5: sdcTable.Rnw:150-151
###################################################
print(lev.V1)


###################################################
### code chunk number 6: sdcTable.Rnw:155-156
###################################################
print(lev.V2)


###################################################
### code chunk number 7: sdcTable.Rnw:160-161
###################################################
print(lev.V3)


###################################################
### code chunk number 8: sdcTable.Rnw:176-177
###################################################
print(tail(completeData))


###################################################
### code chunk number 9: sdcTable.Rnw:180-183
###################################################
levComp.V1 <- as.character(dimList$V1[,2])
levComp.V2 <- as.character(dimList$V2[,2])
levComp.V3 <- as.character(dimList$V3[,2])


###################################################
### code chunk number 10: sdcTable.Rnw:193-194
###################################################
print(levComp.V1)


###################################################
### code chunk number 11: sdcTable.Rnw:198-199
###################################################
print(levComp.V2)


###################################################
### code chunk number 12: sdcTable.Rnw:203-204
###################################################
print(levComp.V3)


###################################################
### code chunk number 13: sdcTable.Rnw:208-210
###################################################
x <- completeData[nrow(completeData),]



###################################################
### code chunk number 14: sdcTable.Rnw:238-240
###################################################
subTots.V1 <- setdiff(levComp.V1, lev.V1)
print(subTots.V1)


###################################################
### code chunk number 15: sdcTable.Rnw:243-245
###################################################
subTots.V2 <- setdiff(levComp.V2, lev.V2)
print(subTots.V2)


###################################################
### code chunk number 16: sdcTable.Rnw:248-250
###################################################
subTots.V3 <- setdiff(levComp.V3, lev.V3)
print(subTots.V3)


###################################################
### code chunk number 17: sdcTable.Rnw:290-293
###################################################
dimV1 <- matrix(nrow=0, ncol=2)
dimV1 <- rbind(dimV1, c('@','Tot'))
print(dimV1)


###################################################
### code chunk number 18: sdcTable.Rnw:306-311
###################################################
mat <- matrix(nrow=4, ncol=2)
mat[,1] <- rep('@@',4)
mat[,2] <- LETTERS[1:4]
dimV1 <- rbind(dimV1, mat)
print(dimV1)


###################################################
### code chunk number 19: sdcTable.Rnw:319-325
###################################################
mat <- matrix(nrow=3, ncol=2)
mat[,1] <- rep('@@@',3)
mat[,2] <- c('Ba','Bb','Bc')

dimV1 <- rbind(dimV1, mat)
print(dimV1)


###################################################
### code chunk number 20: sdcTable.Rnw:335-339
###################################################
dimV1 <- dimV1[c(1:3,6:8, 4:5),]

#dimV1 <- as.data.frame(dimV1, stringsAsFactors=FALSE)
print(dimV1, row.names=FALSE)


###################################################
### code chunk number 21: sdcTable.Rnw:356-362
###################################################
dimV2 <- matrix(nrow=3, ncol=2)
dimV2[,1] <- c('@','@@','@@')
dimV2[,2] <- c('Tot','m','w')
#dimV2 <- mat
#dimV2 <- as.data.frame(mat, stringsAsFactors=FALSE)
print(dimV2, row.names=FALSE)


###################################################
### code chunk number 22: sdcTable.Rnw:377-382
###################################################
dimV3 <- matrix(nrow=7, ncol=2)
dimV3[,1] <- c('@',rep('@@',6))
dimV3[,2] <- c('Tot',letters[1:6])
#dimV3 <- as.data.frame(mat, stringsAsFactors=FALSE)
print(dimV3, row.names=FALSE)


###################################################
### code chunk number 23: sdcTable.Rnw:404-414
###################################################
dimInfo <- list(V1=dimV1, V2=dimV2, V3=dimV3) 

prob.microDat <- makeProblem(
	data=microData, 
	dimList=dimList, 
	dimVarInd=1:3, 
	freqVarInd=NULL, 
	numVarInd=4:5, 
	weightInd=NULL,
	sampWeightInd=NULL) 


###################################################
### code chunk number 24: sdcTable.Rnw:447-458
###################################################
### problem from complete data ###
dimInfo <- list(V1=dimV1, V2=dimV2, V3=dimV3) 
prob.completeDat <- makeProblem(
	data=completeData, 
	dimList=dimList, 
	dimVarInd=1:3, 
	freqVarInd=4, 
	numVarInd=5:6, 
	weightInd=NULL,
	sampWeightInd=NULL) 
#print(table(prob.completeDat@problemInstance@Freq))


###################################################
### code chunk number 25: sdcTable.Rnw:467-468
###################################################
all(c(class(prob.microDat), class(prob.completeDat))=='sdcProblem')


###################################################
### code chunk number 26: sdcTable.Rnw:476-479
###################################################
counts1 <- getInfo(prob.completeDat, type='freq')
counts2 <- getInfo(prob.microDat, type='freq')
all(counts1==counts2)


###################################################
### code chunk number 27: sdcTable.Rnw:492-493
###################################################
length(which(counts1 <= 10) )


###################################################
### code chunk number 28: sdcTable.Rnw:510-511
###################################################
prob.completeDat <- primarySuppression(prob.completeDat, type='freq', maxN=10)


###################################################
### code chunk number 29: sdcTable.Rnw:524-525
###################################################
print(table(getInfo(prob.completeDat, type='sdcStatus')))


###################################################
### code chunk number 30: sdcTable.Rnw:528-529
###################################################
nrPrimSupps <- length(which(getInfo(prob.completeDat, type='sdcStatus')=='u'))


###################################################
### code chunk number 31: sdcTable.Rnw:581-583
###################################################
finalData <- getInfo(resOPT, type='finalData')
print(head(finalData))


###################################################
### code chunk number 32: sdcTable.Rnw:591-592
###################################################
summary(resOPT)


