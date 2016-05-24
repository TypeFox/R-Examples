### R code from vignette source 'other-Bessels.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=75)
library(Bessel)


###################################################
### code chunk number 2: gsl-do
###################################################
library(gsl)


###################################################
### code chunk number 3: gsl-help (eval = FALSE)
###################################################
## library(gsl)
## ?bessel_Knu


###################################################
### code chunk number 4: bessel_Inu_scaled
###################################################
   x <- (1:500)*50000; b2 <- BesselI(x, pi, expo=TRUE)
   b1 <- bessel_Inu_scaled(pi, x)
   all.equal(b1,b2,tol=0) ## "Mean relative difference: 1.544395e-12"

   ## the accuracy is *as* limited (probably):
   b1 <- bessel_Inu_scaled(pi, x, give=TRUE)
   summary(b1$err)


###################################################
### code chunk number 5: bessel_Inu-relErr
###################################################
    range(b1$err/ b1$val)


###################################################
### code chunk number 6: sessionInfo
###################################################
toLatex(sessionInfo())


