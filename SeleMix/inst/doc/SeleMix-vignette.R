### R code from vignette source 'SeleMix-vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: SeleMix-vignette.Rnw:276-278
###################################################
library(Ecdat)
data (Labour)


###################################################
### code chunk number 2: SeleMix-vignette.Rnw:298-299
###################################################
 pairs (log(Labour))


###################################################
### code chunk number 3: SeleMix-vignette.Rnw:307-311
###################################################
options(useFancyQuotes="UTF-8")
#options(useFancyQuotes=FALSE)
options(width=65)
options(warn=-1)


###################################################
### code chunk number 4: SeleMix-vignette.Rnw:333-337
###################################################
library (SeleMix)                          # Load the library
y.names<- c("capital", "output" )          # vector of y variables

est1 <- ml.est(y=Labour[,y.names])         # model estimation


###################################################
### code chunk number 5: SeleMix-vignette.Rnw:342-343
###################################################
head(est1$ypred)     # predicted values


###################################################
### code chunk number 6: SeleMix-vignette.Rnw:346-347
###################################################
head(est1$tau,5)       # vector of posterior probability to be contaminated


###################################################
### code chunk number 7: SeleMix-vignette.Rnw:352-355
###################################################
head(est1$outlier)       # vector of flag: 1=outlier
n.outlier <-sum(est1$outlier)   # numbers of outliers 
n.outlier


###################################################
### code chunk number 8: SeleMix-vignette.Rnw:362-364
###################################################
est1$is.conv       # TRUE convergence is reached
est1$n.iter        # number of iterations


###################################################
### code chunk number 9: SeleMix-vignette.Rnw:369-370
###################################################
est1$bic.aic       # bic and aic 


###################################################
### code chunk number 10: SeleMix-vignette.Rnw:379-383
###################################################
 sel1 <- sel.edit(y=Labour[,y.names], ypred=est1$ypred, t.sel=0.02)
 (n.sel <-sum(sel1[,"sel"]))   # number of influential observations

 head(sel1,3)                # first lines of result matrix


###################################################
### code chunk number 11: SeleMix-vignette.Rnw:387-390
###################################################
  appo <- addmargins( table(outlier=est1$outlier, influential=sel1[,"sel"]) )
  dimnames(appo)  <- list(c("non outliers ", "n.outl", "Sum"),
    c("non influential", "influential errors", "Sum"))


###################################################
### code chunk number 12: SeleMix-vignette.Rnw:393-397
###################################################
 library(xtable)
 xtable(as.matrix(appo), caption="Outliers vs Influential Errors", label="tab:tab1",
        display=rep("d",ncol(appo)+1) )



###################################################
### code chunk number 13: fig1
###################################################
sel.pairs (Labour[,y.names], est1$outlier, sel1[,"sel"])


