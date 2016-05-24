### R code from vignette source 'sss.Rnw'

###################################################
### code chunk number 1: Setup
###################################################
library(sss)
filenameSSS <- "sample.sss"
filenameASC <- "sample.asc"


###################################################
### code chunk number 2: Read
###################################################
sssXML <- readSSSmetadata(filenameSSS)
sss <- parseSSSmetadata(sssXML)
asc <- readSSSdata(filenameASC)


###################################################
### code chunk number 3: Display
###################################################
df <- read.sss(filenameSSS, filenameASC)
print(str(df))
print(df$Q1)


