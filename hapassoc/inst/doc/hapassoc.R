### R code from vignette source 'hapassoc.Rnw'

###################################################
### code chunk number 1: hapassoc.Rnw:411-412 (eval = FALSE)
###################################################
## install.packages('hapassoc')


###################################################
### code chunk number 2: hapassoc.Rnw:418-419
###################################################
library(hapassoc)


###################################################
### code chunk number 3: hapassoc.Rnw:535-536 (eval = FALSE)
###################################################
## help(pre.hapassoc)


###################################################
### code chunk number 4: hapassoc.Rnw:594-595
###################################################
data(hypoDat)


###################################################
### code chunk number 5: hapassoc.Rnw:599-600
###################################################
hypoDat[1:5,]


###################################################
### code chunk number 6: hapassoc.Rnw:619-620
###################################################
example1.haplos <- pre.hapassoc(hypoDat,numSNPs=3)


###################################################
### code chunk number 7: hapassoc.Rnw:625-626
###################################################
example2.haplos <- pre.hapassoc(hypoDat[,c(1:4,7:8)],numSNPs=2)


###################################################
### code chunk number 8: hapassoc.Rnw:630-631
###################################################
example2.haplos <- pre.hapassoc(hypoDat[,1:6],numSNPs=2)


###################################################
### code chunk number 9: hapassoc.Rnw:635-637
###################################################
example2.haplos <- pre.hapassoc(hypoDat[,c(1:2,5:8)],numSNPs=2)
example2.haplos <- pre.hapassoc(hypoDat,numSNPs=2)


###################################################
### code chunk number 10: hapassoc.Rnw:656-657
###################################################
example1.haplos$haploDM[1:7,]


###################################################
### code chunk number 11: hapassoc.Rnw:666-667
###################################################
example1.haplos$nonHaploDM[1:7,]


###################################################
### code chunk number 12: hapassoc.Rnw:675-676
###################################################
example1.haplos$initFreq


###################################################
### code chunk number 13: hapassoc.Rnw:686-688
###################################################
example1.regr1<-hapassoc(affected ~ attr+h000+h010+h011+h100+pooled,
                           example1.haplos,family=binomial())


###################################################
### code chunk number 14: hapassoc.Rnw:692-694 (eval = FALSE)
###################################################
## example1.regr<-hapassoc(affected ~ .,baseline="h001",
##                         example1.haplos,family=binomial())


###################################################
### code chunk number 15: hapassoc.Rnw:708-710
###################################################
example1.regr<-hapassoc(affected == 0 ~ attr+h000+h010+h011+h100+pooled,
                        example1.haplos,family=binomial())


###################################################
### code chunk number 16: hapassoc.Rnw:722-724 (eval = FALSE)
###################################################
## example1.regr2 <- hapassoc(affected ~ attr + I(h000==2) + I(h001==2), 
##                            example1.haplos, family=binomial())


###################################################
### code chunk number 17: hapassoc.Rnw:727-728
###################################################
summary(example1.regr)


###################################################
### code chunk number 18: hapassoc.Rnw:731-733
###################################################
example1.regr2 <- hapassoc(affected ~ attr, example1.haplos, family=binomial())
anova(example1.regr1,example1.regr2)


###################################################
### code chunk number 19: hapassoc.Rnw:741-743
###################################################
data(hypoDatGeno)
hypoDatGeno[1:5,]


###################################################
### code chunk number 20: hapassoc.Rnw:757-761
###################################################
hypoDatGeno[1,"M1"]<-NA #First subject missing genotype at M1
#Add a level 'A' to M2
levels(hypoDatGeno[,"M2"])<-c(levels(hypoDatGeno[,"M2"]),"A")
hypoDatGeno[2,"M2"]<-"A" #Modify second subject's genotype


###################################################
### code chunk number 21: hapassoc.Rnw:768-770
###################################################
example2.haplos<-pre.hapassoc(hypoDatGeno,3,maxMissingGenos=2,allelic=F)
example2.haplos$nonHaploDM[142:148,]


###################################################
### code chunk number 22: hapassoc.Rnw:785-786 (eval = FALSE)
###################################################
## pre.hapassoc(hypoDatGeno,3,maxMissingGenos=0,allelic=F,verbose=FALSE)


###################################################
### code chunk number 23: hapassoc.Rnw:799-800
###################################################
example2.haplos$init


###################################################
### code chunk number 24: hapassoc.Rnw:806-809
###################################################
example2.regr<-hapassoc(attr~hAAA+hACA+hACC+hCAA+pooled,
                        example2.haplos,family=gaussian())
summary(example2.regr)


