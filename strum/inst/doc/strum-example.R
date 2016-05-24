### R code from vignette source 'strum-example.Rnw'

###################################################
### code chunk number 1: strum-example.Rnw:20-21
###################################################
library(strum)


###################################################
### code chunk number 2: strum-example.Rnw:23-24
###################################################
set.seed(1)


###################################################
### code chunk number 3: r2
###################################################
assoForm1 = 
  'L1 =~ P1 + P2 + P3 + <e>
   L1 ~ aSNP + <p,e>
  '
myAssoModel = createStrumModel(formulas = assoForm1)
myAssoModel


###################################################
### code chunk number 4: r3
###################################################
dName = system.file("extdata/example_ped.csv", package = "strum")
dF = read.csv(dName, header=T)[,1:18]
names(dF) = c("family","id", "father","mother",names(dF)[5:18])
myAssoData = createStrumData(dF, "Pedigree")
myAssoData


###################################################
### code chunk number 5: r4
###################################################
myAssoResult = strum(myAssoModel, myAssoData)


###################################################
### code chunk number 6: r5
###################################################
myAssoResult


###################################################
### code chunk number 7: r6
###################################################
linkForm1 = 
  'L1 =~ P1 + P2 + P3 + <e>
   L1 ~ <a,p,e>
  '
myLinkModel = createStrumModel(formulas = linkForm1)
myLinkModel


###################################################
### code chunk number 8: r7
###################################################
iName = system.file("extdata/GENIBD.chr1Ped.ibd", package = "strum")
myLinkData = createStrumData(dF, "Pedigree", ibdFileName=iName)


###################################################
### code chunk number 9: r8 (eval = FALSE)
###################################################
## myLinkResultAll = strum(myLinkModel, myLinkData)


###################################################
### code chunk number 10: r9
###################################################
mNames = c("chr1marker1", "chr1marker2")
myLinkResult = strum(myLinkModel, myLinkData, ibdMarkers=mNames)


###################################################
### code chunk number 11: r10
###################################################
myLinkResult[[1]]


###################################################
### code chunk number 12: r11
###################################################
semForm1 = 
  'bp =~ SBP + DBP
   anger =~ A1 + A2
   stress =~ S1 + S2
   bp ~ anger + stress
   stress ~ anger + rs6040343
   var(stress) = .1
  '
mySemModel = createStrumModel(formulas = semForm1)
mySemModel


###################################################
### code chunk number 13: r12
###################################################
mySemData = createStrumData(dF, "Pedigree")
mySemData


###################################################
### code chunk number 14: r13
###################################################
mySemResult = strum(mySemModel, mySemData)


###################################################
### code chunk number 15: 14
###################################################
mySemResult


###################################################
### code chunk number 16: sessionInfo
###################################################
sessionInfo();


