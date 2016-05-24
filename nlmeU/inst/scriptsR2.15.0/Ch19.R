###################################################
### code chunk: Chap19init
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)

packageVersion("nlmeU")
packageVersion("nlme")
packageVersion("lattice")
packageVersion("reshape")
packageVersion("plyr")

sessionInfo()
  
library(lattice)
data(fcat, package = "nlmeU")

## ---->>>> NOTE: Code used in Panels R19.1 - R19.7 is stored in Ch19mer.R file


###################################################
### code chunk: R19.8
###################################################
library(nlme)  

fcat1 <- within(fcat, one1 <- one2 <- 1L)
system.time(
   fm19.2 <-                                     
      lme(scorec ~ 1, 
          random = list(                      
          one1 = pdIdent(~target - 1), 
          one2 = pdIdent(~id - 1)),
          data = fcat1))
fm19.2                                          # M19.2: (19.2)


###################################################
### code chunk: R19.9a
###################################################
fm19.2$call$data            # Data name
logLik(fm19.2)              # REML value
fixef(fm19.2)               # beta
fm19.2$dims$N               # Number of observations


###################################################
### code chunk: R19.9b
###################################################

## getVarCov(fm19.2)         # Commented out: Not implemented for multiple levels of nesting
VarCorr(fm19.2)             





###################################################
### code chunk: R19.10
###################################################
intervals(fm19.2)


###################################################
### code chunk: R19.11a
###################################################

rnf <- ranef(fm19.2)
## plot(rnf)                  # Commented out: Error in eval(expr, envir, enclos) : object '.pars' not found

###################################################
### code chunk: R19.11b
###################################################
rnft <- lapply(rnf, t)        # Transpose components
names(plxLis  <-              # Auxiliary list ...   
   lapply(rnft, qqnorm,       # ... with two components
          plot.it = FALSE))                
plx <- 
   lapply(plxLis, 
          FUN = function(el) xyplot(y ~ x, data = el, grid = TRUE))
plx[["one1"]]                 # Q-Q plot for id (see Fig. 19.1a)
plx[["one2"]]                 # Q-Q plot for target (see Fig. 19.1b)


###################################################
### code chunk: R19.11c
###################################################
rsd2 <-                               # Equivalent to raw residuals
   resid(fm19.2, type = "pearson") 

  xyplot(rsd2 ~ target, data = fcat1)     # Fig. not shown
   bwplot(rsd2 ~ target, data = fcat1     # Fig. 19.5 
    , # panel = panel.bwxplot             # User defined panel needed (not shown)
)


###################################################
### code chunk: R19.12
###################################################
nItms <- c(4, 6, 8, 5, 9, 6, 8, 6, 5)  # Number of items per target
(nms <- levels(fcat1$target))          # Names extracted ...
names(nItms) <- nms                      # ... and assigned
fcatm <-                               # Add to the data frame...
   within(fcat1, 
          {               
           nItems <- nItms[as.numeric(target)] #... no. of items...
           scorem <- scorec/nItems          # ... mean target-score.
          }) 
(varWghts <- 1/sqrt(nItms))            # Variance function weights
(fxdW <- varWghts[-1]/0.5)             # Ratios wrt the 1st element

fm19.3 <-                              # M19.3
   lme(scorem ~ 1,
   random = list(one1 = pdIdent(~target - 1),
                 one2 = pdIdent(~id - 1)),
   weight = varIdent(form = ~1|target, fixed = fxdW),
   data = fcatm)

###################################################
### code chunk: R19.13
###################################################
summary(fm19.3)$tTable
VarCorr(fm19.3)

sessionInfo()
