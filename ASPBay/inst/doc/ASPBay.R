### R code from vignette source 'ASPBay.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: prompton
###################################################
options(prompt="> ", continue = "+ ");


###################################################
### code chunk number 2: promptoff
###################################################
options(prompt=" ", continue=" ");


###################################################
### code chunk number 3: ASPBay.Rnw:32-33
###################################################
options(prompt="> ", continue = "+ ");


###################################################
### code chunk number 4: desc
###################################################
require(ASPBay)
options(width = 90)
desc <- packageDescription("ASPBay")


###################################################
### code chunk number 5: ASPBay.Rnw:102-103
###################################################
data(ASPData)


###################################################
### code chunk number 6: ASPBay.Rnw:107-110
###################################################
Score <- ASP.Score(ASPData$Control, ASPData$Index, ASPData$IBD)
Score$Value
Score$Pvalue


###################################################
### code chunk number 7: ASPBay.Rnw:115-117
###################################################
Select <- ASP.Selection(ASPData$Control, ASPData$Index, ASPData$IBD)
Select$SNPnames_subset


###################################################
### code chunk number 8: ASPBay.Rnw:121-123
###################################################
M15 <- ASP.Bayesian(1e7, ASPData$Control, ASPData$Index, ASPData$IBD, 15, thin = 10, sd.psi=0.03)
G15 <- Graphs.Bayesian(M15, burn = 1000, print=FALSE)


###################################################
### code chunk number 9: Hexfafb
###################################################
print(G15$hex_fa_fb)


###################################################
### code chunk number 10: Hexr2OR
###################################################
print(G15$hex_r2_OR)


###################################################
### code chunk number 11: HistOR
###################################################
OR <- M15$OR*(M15$OR>=1) + 1/M15$OR*(M15$OR<1)
hist(OR[-(1:10000)], freq=FALSE, breaks=1000, main='Histogram of the OR',
     xlab='OR', xlim=c(1,6))
lines(c(2,2), c(0,100) ,type='l', col='red')


###################################################
### code chunk number 12: Histfa
###################################################
M15$fa <- M15$f_ab + M15$f_aB
fa <- M15$fa*(M15$OR>=1) + (1-M15$fa)*(M15$OR<1)
hist(fa[-(1:10000)], freq=FALSE, breaks=1000,
   	 main='Histogram of the alternative allele\n frequency of the causal variant',
   	 xlab='fa')
lines(c(0.0754,0.0754), c(0,100) ,type='l', col='red')


###################################################
### code chunk number 13: Histfb
###################################################
M15$fb <- M15$f_ab + M15$f_Ab
hist(M15$fb[-(1:10000)], freq=FALSE, breaks=1000,
   	 main='Histogram of the alternative allele\n frequency of the observed variant',
   	 xlab='fb')
lines(c(0.059,0.059), c(0,100) ,type='l', col='red')


###################################################
### code chunk number 14: Histr2
###################################################
M15$D <- M15$f_ab*M15$f_AB - M15$f_aB*M15$f_Ab
M15$r2 <- M15$D**2/( M15$fa*(1-M15$fa)*M15$fb*(1-M15$fb) )
hist(M15$r2[-(1:10000)], freq=FALSE, breaks=1000,
     main='Histogram of linkage disequilibrium', xlab=expression(r^2))
lines(c(0.77,0.77), c(0,100) ,type='l', col='red')


