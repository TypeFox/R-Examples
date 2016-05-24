### R code from vignette source 'indicspeciesTutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: indicspeciesTutorial.Rnw:18-19
###################################################
options(width=67)


###################################################
### code chunk number 2: indicspeciesTutorial.Rnw:25-26
###################################################
library(indicspecies)


###################################################
### code chunk number 3: indicspeciesTutorial.Rnw:33-34
###################################################
data(wetland)


###################################################
### code chunk number 4: indicspeciesTutorial.Rnw:40-42
###################################################
groups = c(rep(1, 17), rep(2, 14), rep(3,10))
groups


###################################################
### code chunk number 5: indicspeciesTutorial.Rnw:45-48
###################################################
wetkm = kmeans(wetland, centers=3)
groupskm = wetkm$cluster
groupskm


###################################################
### code chunk number 6: indicspeciesTutorial.Rnw:59-61
###################################################
indval = multipatt(wetland, groups, 
                   control = how(nperm=999)) 


###################################################
### code chunk number 7: indicspeciesTutorial.Rnw:68-69
###################################################
summary(indval) 


###################################################
### code chunk number 8: indicspeciesTutorial.Rnw:77-78
###################################################
summary(indval, indvalcomp=TRUE)


###################################################
### code chunk number 9: indicspeciesTutorial.Rnw:84-85
###################################################
summary(indval, alpha=1)


###################################################
### code chunk number 10: indicspeciesTutorial.Rnw:88-89
###################################################
indval$sign


###################################################
### code chunk number 11: indicspeciesTutorial.Rnw:95-98
###################################################
wetlandpa = as.data.frame(ifelse(wetland>0,1,0))
phi = multipatt(wetlandpa, groups, func = "r", 
                control = how(nperm=999)) 


###################################################
### code chunk number 12: indicspeciesTutorial.Rnw:103-105
###################################################
phi = multipatt(wetlandpa, groups, func = "r.g", 
                control = how(nperm=999)) 


###################################################
### code chunk number 13: indicspeciesTutorial.Rnw:111-112
###################################################
summary(phi)


###################################################
### code chunk number 14: indicspeciesTutorial.Rnw:117-118
###################################################
round(head(phi$str),3)


###################################################
### code chunk number 15: indicspeciesTutorial.Rnw:121-122
###################################################
round(head(indval$str),3)


###################################################
### code chunk number 16: indicspeciesTutorial.Rnw:132-135
###################################################
indvalori = multipatt(wetland, groups, duleg = TRUE, 
                      control = how(nperm=999)) 
summary(indvalori)


###################################################
### code chunk number 17: indicspeciesTutorial.Rnw:140-143
###################################################
indvalrest = multipatt(wetland, groups, max.order = 2, 
                       control = how(nperm=999)) 
summary(indvalrest)


###################################################
### code chunk number 18: indicspeciesTutorial.Rnw:151-154
###################################################
indvalrest = multipatt(wetland, groups, restcomb = c(1,2,3,5,6), 
                       control = how(nperm=999)) 
summary(indvalrest)


###################################################
### code chunk number 19: indicspeciesTutorial.Rnw:157-158
###################################################
indvalrest$sign


###################################################
### code chunk number 20: indicspeciesTutorial.Rnw:166-168
###################################################
prefstat = strassoc(wetland, cluster=groups, func="A.g")
round(head(prefstat),3)


###################################################
### code chunk number 21: indicspeciesTutorial.Rnw:171-174
###################################################
prefstat = strassoc(wetland, cluster=groups, func="A.g", nboot = 199)
round(head(prefstat$lowerCI),3)
round(head(prefstat$upperCI),3)


###################################################
### code chunk number 22: indicspeciesTutorial.Rnw:181-184
###################################################
prefsign = signassoc(wetland, cluster=groups,  alternative = "two.sided", 
                     control = how(nperm=199)) 
head(prefsign)


###################################################
### code chunk number 23: indicspeciesTutorial.Rnw:192-193
###################################################
coverage(wetland, indvalori)


###################################################
### code chunk number 24: indicspeciesTutorial.Rnw:198-199
###################################################
coverage(wetland, indvalori, At = 0.8)


###################################################
### code chunk number 25: indicspeciesTutorial.Rnw:207-213
###################################################
plotcoverage(x=wetland, y=indvalori, group="1", lty=1)
plotcoverage(x=wetland, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=wetland, y=indvalori, group="3", lty=3, col="red", add=TRUE)
legend(x = 0.01, y=20, 
       legend=c("group 1","group 2", "group 3"),
       lty=c(1,2,3), col=c("black","blue","red"), bty="n")


###################################################
### code chunk number 26: indicspeciesTutorial.Rnw:225-227
###################################################
wetcomb = combinespecies(wetland, max.order = 2)$XC
dim(wetcomb)


###################################################
### code chunk number 27: indicspeciesTutorial.Rnw:230-233
###################################################
indvalspcomb = multipatt(wetcomb, groups, duleg = TRUE, 
                       control = how(nperm=999))
summary(indvalspcomb, indvalcomp = TRUE)


###################################################
### code chunk number 28: indicspeciesTutorial.Rnw:239-242
###################################################
B=strassoc(wetland, cluster=groups ,func="B") 
sel=which(B[,2]>0.2) 
sel


###################################################
### code chunk number 29: indicspeciesTutorial.Rnw:245-247
###################################################
sc= indicators(X=wetland[,sel], cluster=groups, group=2, verbose=TRUE, 
               At=0.5, Bt=0.2)


###################################################
### code chunk number 30: indicspeciesTutorial.Rnw:250-251
###################################################
print(sc, sqrtIVt = 0.6)


###################################################
### code chunk number 31: indicspeciesTutorial.Rnw:259-260
###################################################
coverage(sc)


###################################################
### code chunk number 32: indicspeciesTutorial.Rnw:263-264
###################################################
coverage(sc, At=0.8)


###################################################
### code chunk number 33: indicspeciesTutorial.Rnw:270-274
###################################################
plotcoverage(sc)
plotcoverage(sc, max.order=1, add=TRUE, lty=2, col="red")
legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"), 
       lty=c(1,2), col=c("black","red"), bty="n")


###################################################
### code chunk number 34: indicspeciesTutorial.Rnw:281-283
###################################################
sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)


###################################################
### code chunk number 35: indicspeciesTutorial.Rnw:289-290
###################################################
p<-predict(sc2, wetland)


###################################################
### code chunk number 36: indicspeciesTutorial.Rnw:293-294
###################################################
print(data.frame(groups,p))


