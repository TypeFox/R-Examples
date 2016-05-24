### R code from vignette source 'iqLearn.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: iqLearn.Rnw:477-478
###################################################
options(width=70)


###################################################
### code chunk number 2: iqLearn.Rnw:481-482
###################################################
library (iqLearn)


###################################################
### code chunk number 3: iqLearn.Rnw:485-486
###################################################
data (bmiData)


###################################################
### code chunk number 4: iqLearn.Rnw:491-493
###################################################
dim (bmiData)
head (bmiData)


###################################################
### code chunk number 5: iqLearn.Rnw:497-503
###################################################
bmiData$A1[which (bmiData$A1=="MR")] = 1
bmiData$A1[which (bmiData$A1=="CD")] = -1
bmiData$A2[which (bmiData$A2=="MR")] = 1
bmiData$A2[which (bmiData$A2=="CD")] = -1
bmiData$A1 = as.numeric (bmiData$A1)
bmiData$A2 = as.numeric (bmiData$A2)


###################################################
### code chunk number 6: iqLearn.Rnw:507-509
###################################################
y = -100*(bmiData$month12_BMI -
          bmiData$baseline_BMI)/bmiData$baseline_BMI


###################################################
### code chunk number 7: iqLearn.Rnw:552-555
###################################################
fitIQ2 = learnIQ2 (y ~ gender + parent_BMI + month4_BMI + 
  A2*(parent_BMI + month4_BMI), data=bmiData, treatName="A2", 
  intNames=c ("parent_BMI", "month4_BMI"))


###################################################
### code chunk number 8: iqLearn.Rnw:572-574
###################################################
s2vars = bmiData[, c(1,3,5)]
head (s2vars)


###################################################
### code chunk number 9: iqLearn.Rnw:579-580
###################################################
s2ints = c (2,3)


###################################################
### code chunk number 10: iqLearn.Rnw:583-584
###################################################
fitIQ2 = learnIQ2 (H2=s2vars, Y=y, A2=bmiData$A2, s2ints=s2ints)


###################################################
### code chunk number 11: iqLearn.Rnw:588-589
###################################################
summary (fitIQ2)


###################################################
### code chunk number 12: iqLearn.Rnw:595-596
###################################################
plot (fitIQ2)


###################################################
### code chunk number 13: iqLearn.Rnw:605-606
###################################################
fitIQ2$betaHat20


###################################################
### code chunk number 14: iqLearn.Rnw:609-610
###################################################
fitIQ2$betaHat21


###################################################
### code chunk number 15: iqLearn.Rnw:636-640
###################################################
fitIQ1main = learnIQ1main (~ gender + race + parent_BMI + 
  baseline_BMI + A1*(gender + parent_BMI), data=bmiData, 
  treatName="A1", intNames=c ("gender", "parent_BMI"), s2object=fitIQ2)
summary (fitIQ1main);


###################################################
### code chunk number 16: iqLearn.Rnw:649-651
###################################################
s1vars = bmiData[, 1:4]
head (s1vars)


###################################################
### code chunk number 17: iqLearn.Rnw:657-658
###################################################
s1mainInts = c (1,3)


###################################################
### code chunk number 18: iqLearn.Rnw:661-663
###################################################
fitIQ1main = learnIQ1main (object=fitIQ2, H1Main=s1vars,
    A1=bmiData$A1, s1mainInts=s1mainInts)


###################################################
### code chunk number 19: iqLearn.Rnw:671-672
###################################################
plot (fitIQ1main)


###################################################
### code chunk number 20: iqLearn.Rnw:681-682
###################################################
fitIQ1main$alphaHat0


###################################################
### code chunk number 21: iqLearn.Rnw:685-686
###################################################
fitIQ1main$alphaHat1


###################################################
### code chunk number 22: iqLearn.Rnw:702-708
###################################################
fitIQ1cm = learnIQ1cm (~ gender + race + parent_BMI + 
  baseline_BMI + A1*(gender + parent_BMI + baseline_BMI),
  data=bmiData, treatName="A1", intNames=c ("gender", "parent_BMI",
                                    "baseline_BMI"), 
  s2object=fitIQ2);
summary (fitIQ1cm)


###################################################
### code chunk number 23: iqLearn.Rnw:719-720
###################################################
s1cmInts = c (1,3,4)


###################################################
### code chunk number 24: iqLearn.Rnw:723-725
###################################################
fitIQ1cm = learnIQ1cm (object=fitIQ2, H1CMean=s1vars, A1=bmiData$A1,
    s1cmInts=s1cmInts); 


###################################################
### code chunk number 25: iqLearn.Rnw:730-731
###################################################
plot (fitIQ1cm)


###################################################
### code chunk number 26: iqLearn.Rnw:741-742
###################################################
fitIQ1cm$betaHat10


###################################################
### code chunk number 27: iqLearn.Rnw:745-746
###################################################
fitIQ1cm$betaHat11


###################################################
### code chunk number 28: iqLearn.Rnw:754-755
###################################################
fitIQ1var = learnIQ1var (fitIQ1cm)


###################################################
### code chunk number 29: iqLearn.Rnw:760-761
###################################################
fitIQ1var = learnIQ1var (object=fitIQ1cm, method="homo")


###################################################
### code chunk number 30: iqLearn.Rnw:779-786
###################################################
fitIQ1var = learnIQ1var (~ gender + race + parent_BMI + 
  baseline_BMI + A1*(parent_BMI), data=bmiData, treatName="A1",
    intNames=c ("parent_BMI"), method="hetero", cmObject=fitIQ1cm)  
s1varInts = c (3, 4)
fitIQ1var = learnIQ1var (object=fitIQ1cm, H1CVar=s1vars,
    s1sInts=s1varInts, method="hetero")
summary (fitIQ1var)


###################################################
### code chunk number 31: iqLearn.Rnw:791-792
###################################################
plot (fitIQ1var)


###################################################
### code chunk number 32: iqLearn.Rnw:800-801
###################################################
fitIQ1var$gammaHat0


###################################################
### code chunk number 33: iqLearn.Rnw:804-805
###################################################
fitIQ1var$gammaHat1


###################################################
### code chunk number 34: iqLearn.Rnw:822-823
###################################################
fitResids = iqResids (fitIQ1var)


###################################################
### code chunk number 35: iqLearn.Rnw:827-828
###################################################
plot (fitResids)


###################################################
### code chunk number 36: iqLearn.Rnw:846-848
###################################################
fitIQ1 = learnIQ1 (mainObj=fitIQ1main, cmObj=fitIQ1cm, sigObj=fitIQ1var, 
    dens="nonpar")  


###################################################
### code chunk number 37: iqLearn.Rnw:873-880
###################################################
h1 = c (1, 1, 30, 35)
h1main = h1
h1cm = h1
h1var = h1
optIQ1 = IQ1 (mainObj=fitIQ1main, cmObj=fitIQ1cm, sigObj=fitIQ1var,
    dens="nonpar", h1main=h1main, h1cm=h1cm, h1sig=h1var)
optIQ1


###################################################
### code chunk number 38: iqLearn.Rnw:891-894
###################################################
h2 = c (1, 30, 45);
optIQ2 = IQ2 (fitIQ2, h2);
optIQ2


###################################################
### code chunk number 39: iqLearn.Rnw:911-915
###################################################
fitQ2 = qLearnS2 (H2=s2vars, Y=y, A2=bmiData$A2, s2ints=s2ints);
fitQ2 = qLearnS2 (y ~ gender + parent_BMI + month4_BMI +
  A2*(parent_BMI + month4_BMI), data=bmiData, treatName="A2", 
  intNames=c("parent_BMI", "month4_BMI")); 


###################################################
### code chunk number 40: iqLearn.Rnw:924-929
###################################################
fitQ1 = qLearnS1 (object=fitQ2, H1q=s1vars, A1=bmiData$A1,
    s1ints=c(3,4)); 
fitQ1 = qLearnS1 (~ gender + race + parent_BMI + baseline_BMI +
  A1*(gender + parent_BMI), data=bmiData, treatName="A1", 
    intNames=c ("gender", "parent_BMI"), qS2object=fitQ2);  


###################################################
### code chunk number 41: iqLearn.Rnw:937-941
###################################################
fitQ2$betaHat20
fitQ2$betaHat21
fitQ1$betaHat10
fitQ1$betaHat11


###################################################
### code chunk number 42: iqLearn.Rnw:961-969
###################################################
summary (fitQ1)
h1q = c (1, 1, 30, 35);
optQ1 = qLearnQ1 (fitQ1, h1q);
optQ1
summary (fitQ2)
h2q = c (1, 30, 45);
optQ2 = qLearnQ2 (fitQ2, h2q);
optQ2


###################################################
### code chunk number 43: iqLearn.Rnw:1007-1010
###################################################
estVal = value (d1=fitIQ1$optA1, d2=fitIQ2$optA2, Y=y, A1=bmiData$A1, 
  A2=bmiData$A2)
estVal


