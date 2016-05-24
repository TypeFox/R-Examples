### R code from vignette source 'PRISMA.Rnw'

###################################################
### code chunk number 1: PRISMA.Rnw:14-15
###################################################
if (!exists("PRISMA",.GlobalEnv)) library(PRISMA)  


###################################################
### code chunk number 2: PRISMA.Rnw:78-80
###################################################
data(asap)
asap


###################################################
### code chunk number 3: PRISMA.Rnw:84-85
###################################################
asap$data


###################################################
### code chunk number 4: PRISMA.Rnw:94-95
###################################################
asap$group


###################################################
### code chunk number 5: PRISMA.Rnw:102-104
###################################################
dim(getDuplicateData(asap))
dim(asap$unprocessed)


###################################################
### code chunk number 6: PRISMA.Rnw:110-112
###################################################
asap$duplicatecount
sum(asap$duplicatecount)


###################################################
### code chunk number 7: PRISMA.Rnw:127-129
###################################################
asapDim = estimateDimension(asap)
asapDim


