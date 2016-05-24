### R code from vignette source 'MedRegExample.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: MedRegExample.Rnw:20-21
###################################################
options(width=67)


###################################################
### code chunk number 2: MedRegExample.Rnw:25-26
###################################################
library(vegclust)


###################################################
### code chunk number 3: MedRegExample.Rnw:31-34
###################################################
data(medreg)
class(medreg)
length(medreg)


###################################################
### code chunk number 4: MedRegExample.Rnw:37-38
###################################################
strataUp = c(20,50,100,300,600,1200,2400)


###################################################
### code chunk number 5: MedRegExample.Rnw:41-42
###################################################
strataWidths = c(20,30,50,200,300,600,1200)


###################################################
### code chunk number 6: MedRegExample.Rnw:45-46
###################################################
medreg[[1]]


###################################################
### code chunk number 7: MedRegExample.Rnw:51-52
###################################################
medreg.CAP <- CAP(medreg)


###################################################
### code chunk number 8: MedRegExample.Rnw:55-57
###################################################
class(medreg.CAP)
length(medreg.CAP)


###################################################
### code chunk number 9: MedRegExample.Rnw:60-61
###################################################
medreg.CAP[[1]]


###################################################
### code chunk number 10: MedRegExample.Rnw:65-70
###################################################
plot(medreg.CAP, plots="1", sizes=strataUp, xlab="Height (cm)", 
     ylab="Cumulative percent cover")
legend("topright", col=1:5, lty=1, 
       legend=c("Pines","Oaks","Tall shrubs","Scrubs","Grass"), 
       bty="n")


###################################################
### code chunk number 11: MedRegExample.Rnw:77-79
###################################################
medreg.D = vegdiststruct(medreg.CAP, method="bray", 
                         classWidths=strataWidths)


###################################################
### code chunk number 12: MedRegExample.Rnw:82-83
###################################################
as.matrix(medreg.D)[1,2]


###################################################
### code chunk number 13: MedRegExample.Rnw:86-88
###################################################
medreg.Dsqrt = vegdiststruct(medreg.CAP, method="bray", 
                         classWidths=strataWidths, transform="sqrt")


###################################################
### code chunk number 14: MedRegExample.Rnw:92-99
###################################################
par(mfrow=c(2,1), mar=c(4,5,2,1))
X<-cmdscale(medreg.D, k=2)
plot(X, xlab="MDS 1", ylab="MDS 2", asp=1,
     main="Cover untransformed", cex=0.5)
Xsqrt<-cmdscale(medreg.Dsqrt, k=2)
plot(Xsqrt, xlab="MDS 1", ylab="MDS 2", asp=1,
     main="Cover sqrt-transformed", cex=0.5)


###################################################
### code chunk number 15: MedRegExample.Rnw:106-108
###################################################
nclusters = 6
dnoise = 0.40


###################################################
### code chunk number 16: MedRegExample.Rnw:111-113
###################################################
vc<-vegclustdist(medreg.Dsqrt, mobileMemb = nclusters, 
                 method="HNCdd", dnoise=dnoise, nstart=100)


###################################################
### code chunk number 17: MedRegExample.Rnw:116-118
###################################################
medoids<-vc$mobileCenters
print(medoids)


###################################################
### code chunk number 18: MedRegExample.Rnw:121-123
###################################################
cluster<-defuzzify(vc)$cluster
table(cluster)


###################################################
### code chunk number 19: MedRegExample.Rnw:127-132
###################################################
clNum = as.numeric(as.factor(cluster))
plot(Xsqrt, xlab="MDS 1", ylab="MDS 2", 
     pch=clNum, col=clNum)
legend("topleft", col=1:(nclusters+1), pch=1:(nclusters+1),
       legend=levels(as.factor(cluster)), bty="n")


###################################################
### code chunk number 20: MedRegExample.Rnw:138-140
###################################################
CAPm = CAPmeans(medreg.CAP, as.factor(cluster))
names(CAPm)


###################################################
### code chunk number 21: MedRegExample.Rnw:143-144
###################################################
round(CAPm$M4, dig=1)


###################################################
### code chunk number 22: MedRegExample.Rnw:148-159
###################################################
par(mfrow=c(3,2), mar=c(4,4,3,0))
plot(CAPm, plots="M1", sizes = strataWidths, 
     ylab="Percent cover", main="M1")
plot(CAPm, plots="M2", sizes = strataWidths, main="M2")
plot(CAPm, plots="M3", sizes = strataWidths,  
     ylab="Percent cover", main="M3")
plot(CAPm, plots="M4", sizes = strataWidths, main="M4")
plot(CAPm, plots="M5", sizes = strataWidths, 
     xlab="Height (cm)", ylab="Percent cover", main="M5")
plot(CAPm, plots="M6", sizes = strataWidths, 
     xlab="Height (cm)",  main="M6")


