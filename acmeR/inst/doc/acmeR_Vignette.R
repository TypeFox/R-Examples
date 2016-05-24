## ----echo=FALSE,warning=FALSE,results = "hide"---------------------------
#library(devtools)
#load_all(pkg = '/Users/Jake/Documents/Summer_2015/acme/acmeR')
#library(knitr)
#opts_knit$set(root.dir = '/Users/Jake/Documents/Summer_2015',collapse = TRUE, comment = "#>")
library(acmeR)

## ----eval=FALSE----------------------------------------------------------
#  acme.summary()

## ------------------------------------------------------------------------
## Error: Missing csv file 'fname'.  Example 'altamont.csv' created in
##      working directory.  Try acme.summary('altamont.csv').

## ----eval=FALSE----------------------------------------------------------
#  acme.val <- acme.summary(fname = 'altamont.csv',spec=c("BHCO", "HOWR"))

## ------------------------------------------------------------------------
acme.val <- acme.summary(fname=altamont,spec=c("BHCO","HOWR"))
acme.val

## ---- fig.width=9, fig.height=9------------------------------------------
#Carcass count of 5, with output from acme.summary as parameters
acme.post(C=5, Rstar=acme.val$Rstar, T=acme.val$T, I = acme.val$I, 
        xi = 1/2, lam = 0)

## ---- echo=FALSE---------------------------------------------------------
res <-acme.post(C=5, Rstar=acme.val$Rstar, T=acme.val$T, I = acme.val$I, 
        xi = 1/2, lam = 0, plotit=FALSE)

## ---- fig.height=9,fig.width=9-------------------------------------------
#Carcass count at 3, different hpd credible intervals values
acme.post(C=3,Rstar=acme.val$Rstar, T=acme.val$T, I = acme.val$I, 
        xi = 1/2, lam = 0, gam = c(.9,.95))

## ------------------------------------------------------------------------
acme.table(C=0:5, Rstar=acme.val$Rstar, T=acme.val$T, I = acme.val$I, 
        xi = 1/2, lam = 0, gam=c(.5,.9))

