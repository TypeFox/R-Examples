### R code from vignette source 'strum-intro.Rnw'

###################################################
### code chunk number 1: strum-intro.Rnw:29-30
###################################################
library(strum)


###################################################
### code chunk number 2: strum-intro.Rnw:32-33
###################################################
set.seed(1)


###################################################
### code chunk number 3: r2
###################################################
formulas = 
  'L1 =~ P1 + P2 + P3 + <e>
   L1 ~ aSNP + <p,e>
  '
myModel = createStrumModel(formulas = formulas)
myModel


###################################################
### code chunk number 4: r3
###################################################
inPed = system.file("extdata/example_ped.csv", package = "strum")
dfPed = read.csv(inPed, header=T)[,c(1:6,8:10,17)]
names(dfPed)[1:4] = c("family","id", "father","mother")
myPedData = createStrumData(dfPed, "Pedigree")
myPedData


###################################################
### code chunk number 5: r4
###################################################
iName = system.file("extdata/GENIBD.chr1Ped.ibd", package = "strum")
myPedDataIBD = createStrumData(dfPed, "Pedigree", ibdFileName=iName)
myPedDataIBD


###################################################
### code chunk number 6: r5
###################################################
inRaw = system.file("extdata/example_raw.csv", package = "strum")
dfRaw = read.csv(inRaw, header=T)
head(dfRaw)
myRawData = createStrumData(dfRaw, "RawData")
myRawData


###################################################
### code chunk number 7: r6
###################################################
myFitResult = strum(myModel, myPedData)


###################################################
### code chunk number 8: r7 (eval = FALSE)
###################################################
## mNames = c("chr1marker1", "chr1marker2")
## myLinkResult = strum(myLinkModel, myPedIBD, ibdMarkers=mNames)


###################################################
### code chunk number 9: r8 (eval = FALSE)
###################################################
## hap20 = importHapmapData(20)


###################################################
### code chunk number 10: r9
###################################################
#hap20snp10 = hap20[(1:10)*10,]
#save(hap20snp10,file="hap20snp10.Rdata")
#using locally saved copy
inHap = system.file("extdata/hap20snp10.Rdata", package = "strum")
load(file=inHap)
snpStrumMarker = createStrumMarker(hapMapData=hap20snp10)
#snpStrumMarker


###################################################
### code chunk number 11: r10
###################################################
simform = 
  'L1 =~ X1 + 2*X2 + 0.5*X3 + <e>
   L1 ~ aSNP + <p,e>
  '
mySimModel = createSimModel(formulas = simform,
                            tMissingRate = c(0.1),
                            markerInfo = snpStrumMarker)
#mySimModel


###################################################
### code chunk number 12: r11
###################################################
mySimData = simulateStrumData(mySimModel, myPedData)
#mySimData


###################################################
### code chunk number 13: r12
###################################################
simform1 = 'z1 =~ X1 + 0.8*X2 + 0.5*X3 + y'
mySimModel1 = createSimModel(formulas = simform1,
                             defaultError='<e>')
#mySimModel1
mySimData1 = simulateStrumData(mySimModel1, N=150)
#mySimData1


###################################################
### code chunk number 14: r13 (eval = FALSE)
###################################################
## testform = 'z1 =~ X1 + X2 + X3 + y'
## myTestModel = createStrumModel(formulas = testform, defaultError='<e>')
## mySimResult = strum(myTestModel, mySimData1)
## #mySimResult


###################################################
### code chunk number 15: r14
###################################################
myAStrumModel = createStrumModel(formulas = formulas,
                                 ascertainment="proband")
myAStrumModel


###################################################
### code chunk number 16: r15
###################################################
aFunction = function(thisFam)
            {
              aff = (thisFam$disease == 1)
              ascertained = any(aff)
              proband = rep(FALSE, nrow(thisFam))
              if(ascertained)
                pPos = which.min(thisFam$disease == 1)
                proband[pPos] = TRUE
              return(list(aStatus=ascertained, pStatus=proband))
            }

myASimModel = createSimModel(formulas = simform,
                             markerInfo = snpStrumMarker,
                             ascertainment = aFunction)
#myASimModel
myASimData = simulateStrumData(myASimModel, myPedData)
#myASimData


###################################################
### code chunk number 17: sessionInfo
###################################################
sessionInfo();


