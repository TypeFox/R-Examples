### R code from vignette source 'Guide.Stex'

###################################################
### code chunk number 1: Guide.Stex:6-7
###################################################
 options(continue="  ")


###################################################
### code chunk number 2: Guide.Stex:23-24
###################################################
library("tsfa") 


###################################################
### code chunk number 3: Guide.Stex:36-62
###################################################

data("CanadianMoneyData.asof.28Jan2005", package="CDNmoney")
data("CanadianCreditData.asof.28Jan2005", package="CDNmoney")

cpi <- 100 * M1total / M1real
seriesNames(cpi) <- "CPI"
popm <- M1total / M1PerCapita
seriesNames(popm) <- "Population of Canada"

z <- tframed(tbind(
    MB2001,
    MB486 + MB452 + MB453 ,
    NonbankCheq,
    MB472 + MB473 + MB487p,
    MB475,
    NonbankNonCheq + MB454 + NonbankTerm + MB2046 + MB2047 + MB2048 +
    MB2057 + MB2058 + MB482),
    names=c("currency", "personal cheq.", "NonbankCheq",
    "N-P demand & notice", "N-P term", "Investment" )
  )


TotalMoney <- tframed(rowSums(z), tframe(z))

z <- tbind (z, ConsumerCredit, ResidentialMortgage,
    ShortTermBusinessCredit, OtherBusinessCredit)


###################################################
### code chunk number 4: Guide.Stex:67-72
###################################################
z <-tfwindow(z, start=c(1981,11), end=c(2004,11))

scale <- tfwindow(1e8 /(popm * cpi), tf=tframe(z))

MBandCredit <- sweep(z, 1, scale, "*")


###################################################
### code chunk number 5: Guide.Stex:79-80
###################################################
tfplot(MBandCredit, graphs.per.page=3)


###################################################
### code chunk number 6: Guide.Stex:83-84
###################################################
tfplot(diff(MBandCredit), graphs.per.page=3)


###################################################
### code chunk number 7: Guide.Stex:89-93
###################################################
start(MBandCredit)
end(MBandCredit)
Tobs(MBandCredit)
DX <- diff(MBandCredit, lag=1)


###################################################
### code chunk number 8: Guide.Stex:96-97
###################################################
Tobs(MBandCredit)


###################################################
### code chunk number 9: Guide.Stex:100-101
###################################################
nseries(MBandCredit)


###################################################
### code chunk number 10: Guide.Stex:104-105
###################################################
colMeans(DX)


###################################################
### code chunk number 11: Guide.Stex:108-109
###################################################
sqrt(diag(cov(DX)))


###################################################
### code chunk number 12: Guide.Stex:116-118
###################################################
zz <- eigen(cor(diff(MBandCredit, lag=1)), symmetric=TRUE)[["values"]]
print(zz)


###################################################
### code chunk number 13: Guide.Stex:121-123
###################################################
par(omi=c(0.1,0.1,0.1,0.1),mar=c(4.1,4.1,0.6,0.1))
plot(zz, ylab="Value", xlab="Eigenvalue Number", pch=20:20,cex=1,type="o")


###################################################
### code chunk number 14: Guide.Stex:129-132
###################################################
z <- FAfitStats(MBandCredit)    
print(z, digits=3)
c2withML <- estTSF.ML(MBandCredit, 2)


###################################################
### code chunk number 15: Guide.Stex:146-157
###################################################
z <- matrix(0,10,3)
z[matrix(c( 1,6,2,1:3),3,2)] <- c(10, 56, 41)
c3withML <- estTSF.ML(MBandCredit, 3, BpermuteTarget=z)

z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(13, 54, 37, 24)
c4withML  <- estTSF.ML(MBandCredit, 4, BpermuteTarget=z)

z <- matrix(0,10,5)
z[matrix(c( 1,6,2,7,5,1:5),5,2)] <- c(13, 67, 34, 30, 72)
c5withML <- estTSF.ML(MBandCredit, 5, BpermuteTarget=z)


###################################################
### code chunk number 16: Guide.Stex:162-163
###################################################
print(DstandardizedLoadings(c4withML) )


###################################################
### code chunk number 17: Guide.Stex:167-168
###################################################
print(c4withML$Phi, digits=3)


###################################################
### code chunk number 18: Guide.Stex:173-174
###################################################
print(1 - c4withML$stats$uniquenesses)


###################################################
### code chunk number 19: Guide.Stex:179-182
###################################################
print(1 - c2withML$stats$uniquenesses)
print(1 - c3withML$stats$uniquenesses)
print(1 - c5withML$stats$uniquenesses)


###################################################
### code chunk number 20: Guide.Stex:187-188
###################################################
print(loadings(c4withML) )


###################################################
### code chunk number 21: Guide.Stex:193-200
###################################################
tfplot(ytoypc(factors(c4withML)),
       Title= "Factors from 4 factor model (year-to-year growth rate)",
       lty=c("solid"),
       col=c("black"),
       xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 22: Guide.Stex:205-212
###################################################
tfplot(factors(c4withML),
       Title="Factors from 4 factor model",
       lty=c("solid"),
       col=c("black"),
       xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 23: Guide.Stex:219-232
###################################################

z <- explained(c4withML)

tfplot(ytoypc(MBandCredit), ytoypc(z), series=1:5, graphs.per.page=5,
       lty=c("solid", "dashed"), 
       col=c("black", "red"),
       ylab=c("currency", "personal cheq.", "NonbankCheq",
              "N-P demand & notice", "N-P term"),
       ylim=list(NULL,NULL,c(-100,100),NULL,NULL),
       Title=
         "Explained money indicator 1-5 (year-to-year growth rate) using 4 factors",
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 24: Guide.Stex:237-248
###################################################
tfplot(ytoypc(MBandCredit), ytoypc(explained(c4withML)), series=6:10,
       graphs.per.page=5,
       lty=c("solid", "dashed"), 
       col=c("black", "red"),
       ylab=c("","","","","",
              "Investment","Consumer Credit", "Residential Mortgage",
              "Short Term Business Credit", "Other Business Credit"),
       Title=
         "Explained money indicator 6-10 (year-to-year growth rate)using 4 factors",
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 25: Guide.Stex:253-260
###################################################

tfplot(  MBandCredit,  explained(c4withML),  series=1:5, graphs.per.page=5,
        lty=c("solid", "dashed"), 
        col=c("black", "red"),
        Title= "Explained money indicators 1-5 using 4 factors",
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 26: Guide.Stex:265-271
###################################################
tfplot(  MBandCredit,  explained(c4withML),  series=6:10, graphs.per.page=5,
        lty=c("solid", "dashed"), 
       col=c("black", "red"),
        Title= "Explained money indicator 6-10 using 4 factors",
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 27: Guide.Stex:274-275
###################################################
tfplot(  diff(MBandCredit), diff(explained(c4withML)),   graphs.per.page=3)


###################################################
### code chunk number 28: Guide.Stex:278-279
###################################################
summary(MBandCredit)


###################################################
### code chunk number 29: Guide.Stex:287-298
###################################################
DstandardizedLoadings(c2withML)
print(c2withML$Phi, digits=3)

DstandardizedLoadings(c3withML)
print(c3withML$Phi, digits=3)

print(DstandardizedLoadings(c5withML), digits=3)
print(c5withML$Phi, digits=3)

DstandardizedLoadings(c4withML)
print(c2withML$Phi, digits=3)


###################################################
### code chunk number 30: Guide.Stex:302-312
###################################################
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c2withML)),
       ytoypc(factors(c3withML)),
       ytoypc(factors(c5withML)), series=1:2,
       xlab=c(""),ylab=c("factor 1","factor 2"),
       lty=c("solid", "dotdash", "dashed", "dotted"),
       col=c("black","green","red","blue"),
       Title= paste("Factors transaction and long ",
        "(year-to-year growth rate) using 2, 3, 4 and 5 factor models", sep=""),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 31: Guide.Stex:316-326
###################################################
tfplot(ytoypc(factors(c4withML)),
       ytoypc(factors(c3withML)),
       ytoypc(factors(c5withML)), series=3,
       lty=c("solid", "dashed", "dotted"),
       xlab=c(""),ylab=c("","","factor 3"),
       col=c("black","red","blue"),
       Title= paste("Factor near ",
        "(year-to-year growth rate) using 3, 4 and 5 factor models", sep=""),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 32: Guide.Stex:335-340
###################################################
c4withMLg0.5 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
     rotation="oblimin", rotationArgs=list(gam=0.5))
loadings(c4withMLg0.5)
DstandardizedLoadings(c4withMLg0.5)
DstandardizedLoadings(c4withMLg0.5) - DstandardizedLoadings(c4withML)


###################################################
### code chunk number 33: Guide.Stex:343-344
###################################################
summary(c4withMLg0.5)


###################################################
### code chunk number 34: Guide.Stex:349-369
###################################################
c4withMLgneg0.5 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
     rotation="oblimin", rotationArgs=list(gam=-0.5))
loadings(c4withMLgneg0.5)
DstandardizedLoadings(c4withMLgneg0.5)
DstandardizedLoadings(c4withMLgneg0.5) - DstandardizedLoadings(c4withML)
summary(c4withMLgneg0.5)

c4withMLgneg1.0 <- estTSF.ML(MBandCredit, 4, BpermuteTarget=loadings(c4withML),
     rotation="oblimin", rotationArgs=list(gam=-1.0))
loadings(c4withMLgneg1.0)
DstandardizedLoadings(c4withMLgneg1.0)
DstandardizedLoadings(c4withMLgneg1.0) - DstandardizedLoadings(c4withML)
summary(c4withMLgneg1.0)

c4withMLbQ <- estTSF.ML(MBandCredit, 4, rotation="bentlerQ",
     BpermuteTarget=loadings(c4withML))
loadings(c4withMLbQ)
DstandardizedLoadings(c4withMLbQ)
DstandardizedLoadings(c4withMLbQ) - DstandardizedLoadings(c4withML)
summary(c4withMLbQ)


###################################################
### code chunk number 35: Guide.Stex:374-385
###################################################
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLg0.5)),
       ytoypc(factors(c4withMLgneg0.5)), ytoypc(factors(c4withMLgneg1.0)),
       ytoypc(factors(c4withMLbQ)),
       xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
       lty=c("solid", "dashed", "dotted", "dotdash", "longdash"),
       col=c("black","red","blue","green","pink"),
       Title= paste(
         "Factors from various 4 factor models (year-to-year growth rate)",
         "\n and oblimin with gam=0 (solid)"),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 36: Guide.Stex:393-397
###################################################
c4withMLgm <- estTSF.ML(MBandCredit, 4, rotation="geominQ",
     BpermuteTarget=loadings(c4withML))
loadings(c4withMLgm)
DstandardizedLoadings(c4withMLgm)


###################################################
### code chunk number 37: Guide.Stex:400-401
###################################################
DstandardizedLoadings(c4withMLgm) - DstandardizedLoadings(c4withML)


###################################################
### code chunk number 38: Guide.Stex:404-405
###################################################
summary(c4withMLgm)


###################################################
### code chunk number 39: Guide.Stex:409-422
###################################################
tfplot(ytoypc(factors(c4withML)), ytoypc(factors(c4withMLgm)),
       xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
       lty=c("solid", "dashed"),
       col=c("black","red"),
       Title= paste(
     "Factors from geomin (dashed) 4 factor model (year-to-year growth rate)",
     "\n and oblimin with gam=0 (solid)"),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)

c4withMLnotNorm  <- estTSF.ML(MBandCredit, 4, normalize=FALSE,
       BpermuteTarget=loadings(c4withML))



###################################################
### code chunk number 40: Guide.Stex:427-430
###################################################
DstandardizedLoadings(c4withML)
DstandardizedLoadings(c4withMLnotNorm)
DstandardizedLoadings(c4withML) - DstandardizedLoadings(c4withMLnotNorm)


###################################################
### code chunk number 41: Guide.Stex:439-457
###################################################
z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(11, 104, 20, 13)
c4withMLbefore90 <- estTSF.ML(tfwindow(MBandCredit, end=c(1989,12)), 4,
         BpermuteTarget=z)

c4withMLafter95 <- estTSF.ML(tfwindow(MBandCredit, start=c(1995,1)), 4,
         BpermuteTarget=loadings(c4withML))

z <- matrix(0,10,4)
z[matrix(c( 1,6,2,7,1:4),4,2)] <- c(11, 104, 20, 13)
c4withMLbefore95 <- estTSF.ML(tfwindow(MBandCredit, end=c(1994,12)), 4,
         BpermuteTarget=z)

c4withMLafter00 <- estTSF.ML(tfwindow(MBandCredit, start=c(2000,1)), 4,
         BpermuteTarget=loadings(c4withML))

c4withML90to00 <- estTSF.ML(tfwindow(MBandCredit, start=c(1990,1), end=c(2000,1)), 4,
         BpermuteTarget=loadings(c4withML))


###################################################
### code chunk number 42: Guide.Stex:462-476
###################################################
tfplot(ytoypc(factors(c4withML)),         ytoypc(factors(c4withMLbefore90)),
       ytoypc(factors(c4withMLbefore95)), ytoypc(factors(c4withMLafter95)),
       ytoypc(factors(c4withMLafter00)),  ytoypc(factors(c4withML90to00)),
       xlab=c(""),ylab=c("factor 1","factor 2","factor 3","factor 4"),
       ylim=list(NULL,c(-20,20),c(-25,40),NULL),
       graphs.per.page=4,
       lty=c("dashed", "dotted", "dotdash", "longdash", "dotted",
             "twodash"),
       col=c("red","blue","green","pink","violet","brown"),
       Title= paste(
           "Factors (year to year growth) using full sample and sub-samples\n",
           "ML estimation with quartimin rotation objective", sep=""),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 43: Guide.Stex:481-495
###################################################
tfplot(ytoypc(MBandCredit),                ytoypc(explained(c4withML)),
       ytoypc(explained(c4withMLbefore90)), ytoypc(explained(c4withMLbefore95)),
       ytoypc(explained(c4withMLafter95)),  ytoypc(explained(c4withMLafter00)),
       ytoypc(explained(c4withML90to00)), series=1:5, graphs.per.page=5,
       ylab=c("currency", "personal cheq.", "NonbankCheq",
              "N-P demand & notice", "N-P term"),
       ylim=list(NULL,NULL,c(-70,70),NULL,c(-70,70)),
       lty=c("solid", "dashed", "dotted", "dotdash", "longdash", "dotted",
             "twodash"),
       col=c("black", "red","blue","green","pink","violet","brown"),
       Title= paste("Explained money indicators 1-5 (year to year growth)\n",
               "using 4 factors, full sample and sub-samples", sep=""),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 44: Guide.Stex:501-515
###################################################
tfplot(ytoypc(MBandCredit), ytoypc(explained(c4withML)),
       ytoypc(explained(c4withMLbefore90)), ytoypc(explained(c4withMLbefore95)),
       ytoypc(explained(c4withMLafter95)),  ytoypc(explained(c4withMLafter00)),
       ytoypc(explained(c4withML90to00)), series=6:10, graphs.per.page=5,
       ylab=c("","","","","",
              "Investment","Consumer Credit", "Residential Mortgage",
              "Short Term Business Credit", "Other Business Credit"),
       lty=c("solid", "dashed", "dotted", "dotdash", "longdash", "dotted",
             "twodash"),
       col=c("black", "red","blue","green","pink","violet","brown"),
       Title= paste("Explained money indicators 6-10 (year to year growth)\n",
               "using 4 factors, full sample and sub-samples", sep=""),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


###################################################
### code chunk number 45: Guide.Stex:524-539
###################################################
M1   <-  tfwindow(M1total, start=c(1981,11), end=c(2004,11)) * scale
seriesNames(M1) <- "Real per Capita M1"

z <-  tframed(MB2001 + MB486 + MB487p + MB452 + MB452adj + MB472 + NonbankCheq)

M1p   <-  tfwindow(z, start=c(1981,11), end=c(2004,11)) * scale
seriesNames(M1p) <- "Real per Capita M1+"

M2pp <-  tfwindow(M1total
                 + MB472 + MB473 + MB452 + MB453 + MB454
                 + NonbankCheq + NonbankNonCheq + NonbankTerm +
                 + MB2046 + MB2047 + MB2048
                 + MB2057 + MB2058, start=c(1981,11), end=c(2004,11))* scale

seriesNames(M2pp) <- "Real per Capita M2++"


###################################################
### code chunk number 46: Guide.Stex:543-547
###################################################
f <- tframed(factors(c4withML)[,1:2], tf=tframe(factors(c4withML)))
mnF <- colMeans(f)
mnM <- colMeans(cbind(M1p, M2pp))
f <- sweep(f, 2, mnM/mnF, "*")


###################################################
### code chunk number 47: Guide.Stex:553-562
###################################################
tfplot(ytoypc(f), ytoypc(cbind(M1, M2pp)), graphs.per.page=2,
       lty=c("dashed", "solid"),
       col=c("red","black"),
       Title=
   paste("(year to year growth) M1+ and M2++ (solid) and scaled Bartlett Predictors\n",
           "computed using ordinary ML parameters (dashed)", sep=""),
       ylab=c("M1+ vs. factor 1", "M2++ vs. factor 2" ),
       par=list(mar=c(2.1, 4.1, 1.1, 0.1)),
       reset.screen=TRUE)


