###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###

#  Remarks and disclaimers

#  These R commands are either those in this book, or designed to
#  otherwise illustrate how R can be used in the analysis of functional
#  data.
#  We do not claim to reproduce the results in the book exactly by these
#  commands for various reasons, including:
#    -- the analyses used to produce the book may not have been
#       entirely correct, possibly due to coding and accuracy issues
#       in the functions themselves
#    -- we may have changed our minds about how these analyses should be
#       done since, and we want to suggest better ways
#    -- the R language changes with each release of the base system, and
#       certainly the functional data analysis functions change as well
#    -- we might choose to offer new analyses from time to time by
#       augmenting those in the book
#    -- many illustrations in the book were produced using Matlab, which
#       inevitably can imply slightly different results and graphical
#       displays
#    -- we may have changed our minds about variable names.  For example,
#       we now prefer "yearRng" to "yearRng" for the weather data.
#    -- three of us wrote the book, and the person preparing these scripts
#       might not be the person who wrote the text
#  Moreover, we expect to augment and modify these command scripts from time
#  to time as we get new data illustrating new things, add functionality
#  to the package, or just for fun.

###
### ch. 9.  Functional Linear Models for Scalar Response
###

#  load the fda package

library(fda)

#  display the data files associated with the fda package

data(package='fda')

#  start the HTML help system if you are connected to the Internet, in
#  order to open the R-Project documentation index page in order to obtain
#  information about R or the fda package.

help.start()

##
## Section 9.1 Functional Linear regression with a Scalar response
##

#  (no computations in this section)

##
## Section 9.2 A Scalar Response Model for Log Annual Precipitation
##

annualprec   = log10(apply(daily$precav,2,sum))

tempbasis65  = create.fourier.basis(c(0,365),65)
tempSmooth65 = smooth.basis(day.5, daily$tempav, tempbasis65)
tempfd65     = tempSmooth65$fd

##
## Section 9.3 Setting Up the Functional Linear Model
##
#  (no computations in this section)

##
## Section 9.4 Three Estimates of the Regression Coefficient
##             Predicting Annual Precipitation
##

templist      = vector("list",2)
templist[[1]] = rep(1,35)
templist[[2]] = tempfd65

# 9.4.1 Low Dimensional Regression Coefficient Function beta

conbasis   = create.constant.basis(c(0,365))
betabasis5 = create.fourier.basis(c(0,365),5)
betalist1  = vector("list",2)
betalist1[[1]] = conbasis
betalist1[[2]] = betabasis5

fRegressList1 = fRegress(annualprec,templist,betalist1)

betaestlist1  = fRegressList1$betaestlist
tempbetafd1   = betaestlist1[[2]]$fd

# Figure 9.1

plot(tempbetafd1, xlab="Day", ylab="Beta for temperature")

coef(betaestlist1[[1]])
# 0.0095 as in the book

annualprechat1 = fRegressList1$yhatfdobj
annualprecres1 = annualprec - annualprechat1
SSE1.1  = sum(annualprecres1^2)
SSE0    = sum((annualprec - mean(annualprec))^2)
(RSQ1   = (SSE0-SSE1.1)/SSE0)
# 0.80 as in the book
(Fratio1 = ((SSE0-SSE1.1)/5)/(SSE1.1/29))
# 22.6 as in the book

# 9.4.2 Coefficient beta Estimate Using a Roughness Penalty

Lcoef = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))

# refit with 35 terms rather than 5 in the fourier basis

betabasis35 = create.fourier.basis(c(0, 365), 35)
lambda      = 10^12.5
betafdPar.  = fdPar(betabasis35, harmaccelLfd, lambda)

betalist2      = betalist1
betalist2[[2]] = betafdPar.

annPrecTemp    = fRegress(annualprec, templist, betalist2)
betaestlist2   = annPrecTemp$betaestlist
annualprechat2 = annPrecTemp$yhatfdobj

print(annPrecTemp$df)

SSE1.2 = sum((annualprec-annualprechat2)^2)
(RSQ2 = (SSE0 - SSE1.2)/SSE0)
# 0.75 as in the book

(Fratio2 = ((SSE0-SSE1.2)/3.7)/(SSE1.2/30.3))
# 25.1 as in the book

# Figure 9.2

plot(annualprechat2, annualprec, lwd=2)
abline(lm(annualprec~annualprechat2), lty='dashed', lwd=2)

# Figure 9.3
# ... see section 9.4.4 below ...

plot(betaestlist2[[2]]$fd, lwd=2)

# Compare with the constant fit:

betalist      = betalist1
betalist[[2]] = fdPar(conbasis)
fRegressList  = fRegress(annualprec, templist, betalist)
betaestlist   = fRegressList$betaestlist

annualprechat = fRegressList$yhatfdobj
SSE1 = sum((annualprec-annualprechat)^2)

(RSQ = (SSE0 - SSE1)/SSE0)
# 0.49 as in the book

(Fratio = ((SSE0-SSE1)/1)/(SSE1/33))
# 31.3 as in the book

# 9.4.3 Choosing Smoothing Parameters

loglam = seq(5,15,0.5)
nlam   = length(loglam)
SSE.CV = rep(NA,nlam)
for (ilam in 1:nlam) {
  print(paste("log lambda =", loglam[ilam]))
  lambda     = 10^(loglam[ilam])
  betalisti  = betalist2
  betafdPar2 = betalisti[[2]]
  betafdPar2$lambda = lambda
  betalisti[[2]] = betafdPar2
  fRegi          = fRegress.CV(annualprec, templist, betalisti)
  SSE.CV[ilam]   = fRegi$SSE.CV
}

#  Figure 9.4

plot(loglam, SSE.CV, type="b", lwd=2,
     xlab="log smoothing parameter lambda",
     ylab="Cross-validation score", cex=1.2)

# 9.4.4 Confidence Intervals

resid   = annualprec - annualprechat2
SigmaE. = sum(resid^2)/(35-annPrecTemp$df)
SigmaE  = SigmaE.*diag(rep(1,35))
y2cMap  = tempSmooth65$y2cMap

stderrList = fRegress.stderr(annPrecTemp, y2cMap, SigmaE)

betafdPar      = betaestlist2[[2]]
betafd         = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd   = betastderrList[[2]]

# Figure 9.3

plot(betafd, xlab="Day", ylab="Temperature Reg. Coeff.",
     ylim=c(-6e-4,1.2e-03), lwd=2)
lines(betafd+2*betastderrfd, lty=2, lwd=1)
lines(betafd-2*betastderrfd, lty=2, lwd=1)

# Section 9.4.5 Scalar Response Models by Functional Principal Components

daybasis365   = create.fourier.basis(c(0, 365), 365)
lambda        = 1e6
tempfdPar365  = fdPar(daybasis365, harmaccelLfd, lambda)
tempSmooth365 = smooth.basis(day.5, daily$tempav,
                              tempfdPar365)
tempfd = tempSmooth365$fd

lambda    = 1e0
tempfdPar = fdPar(daybasis365, harmaccelLfd, lambda)
temppca   = pca.fd(tempfd, 4, tempfdPar)
harmonics = temppca$harmonics

pcamodel = lm(annualprec~temppca$scores)
pcacoefs = summary(pcamodel)$coef
betafd   = pcacoefs[2,1]*harmonics[1] + pcacoefs[3,1]*harmonics[2] +
           pcacoefs[4,1]*harmonics[3]
coefvar  = pcacoefs[,2]^2
betavar  = coefvar[2]*harmonics[1]^2 + coefvar[3]*harmonics[2]^2 +
           coefvar[4]*harmonics[3]^2

# Figure 9.5

plot(betafd, xlab="Day", ylab="Regression Coef.",
     ylim=c(-6e-4,1.2e-03), lwd=2)
s.betavar.2 <- 2*sqrt(betavar)
lines(betafd+s.betavar.2, lty=2, lwd=1)
lines(betafd-s.betavar.2, lty=2, lwd=1)

##
## Section 9.5 Statistical Tests
##

F.res = Fperm.fd(annualprec, templist, betalist)

F.res$Fobs
F.res$qval

##
## Section 9.6 Some Things to Try
##
# (exercises for the reader)

##
## Section 9.7  More to Read
##
