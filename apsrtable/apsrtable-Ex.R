pkgname <- "apsrtable"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('apsrtable')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("apsrtable")
### * apsrtable

flush(stderr()); flush(stdout())

### Name: apsrtable
### Title: APSR-style latex tables with multiple models
### Aliases: apsrtable

### ** Examples
 
     ## Use the example from lm() to show both models:
     ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
     ## Page 9: Plant Weight Data.
     ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
     trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
     group <- gl(2,10,20, labels=c("Ctl","Trt"))
     weight <- c(ctl, trt)
     lm.D9 <- lm(weight ~ group)
     glm.D9 <- glm(weight~group)
     lm.D90 <- lm(weight ~ group - 1) # omitting intercept
     apsrtable(lm.D90, lm.D9, glm.D9, digits=1, align="center", 
               stars="default", model.counter=0, order="rl")
     ## Not run: 
##D apsrtable(lm.D90, lm.D9, glm.D9, digits=1, align="l", 
##D           stars=1, model.counter=0, order="rl",
##D           coef.rows=1, col.hspace="3em", float="sidewaystable")
##D 
##D ## Omit rows by regular expressions
##D apsrtable(lm.D9, omitcoef=expression(grep("\\(",coefnames)))
##D apsrtable(lm.D90,lm.D9,
##D           omitcoef=list("groupCtl",
##D             expression(grep("\\(",coefnames,value=TRUE))
##D             )
##D           )
## End(Not run)




cleanEx()
nameEx("apsrtableSummary")
### * apsrtableSummary

flush(stderr()); flush(stdout())

### Name: apsrtableSummary
### Title: Custom summary functions for output tables
### Aliases: apsrtableSummary apsrtableSummary.gee apsrtableSummary.lrm

### ** Examples

### summary.gee produces z scores but not Pr(z). This converts the relevant columns
### to Pr(z) so that apsrstars() works on it, and places the vector of robust se's in 
### an $se position which apsrtable expects.

apsrtableSummary.gee <- function(x) {
  s <- summary(x)
  newCoef <- coef(s)
  ## which columns have z scores? (two of them in robust case)
  zcols <- grep("z",colnames(newCoef))
  newCoef[,zcols] <- pnorm(abs(newCoef[,zcols]), lower.tail=FALSE)
  colnames(newCoef)[zcols] <- "Pr(z)"
  s$coefficients <- newCoef
  ## put the robust se in $se so that notefunction works automatically
  ## the se checker will overwrite [,4] with pt, but this doesn't matter
  ## because the last column Pr(z) is used by apsrstars() anyway
  ## and the se are pulled from $se.
  if( class(x) == "gee.robust") {
    s$se <- coef(s)[,4]
  }
  return(s)
}



cleanEx()
nameEx("modelInfo")
### * modelInfo

flush(stderr()); flush(stdout())

### Name: modelInfo
### Title: Model fit and diagnostic functions for output
### Aliases: modelInfo modelInfo,summary.lm-method
###   modelInfo,summary.glm-method modelInfo,summary.svyglm-method
###   modelInfo,summary.tobit-method modelInfo,summary.gee-method
###   modelInfo,summary.coxph-method modelInfo,summary.clogit-method
###   modelInfo,summary.negbin-method modelInfo,summary.lrm-method

### ** Examples
 

setMethod("modelInfo", "summary.lm", function(x) {
  env <- sys.parent()
  digits <- evalq(digits, env)
  model.info <- list(
                     "$N$"=formatC(sum(x$df[1:2]),format="d"),
                     "Resid. sd" = formatC(x$sigma,format="f",digits=digits))
  class(model.info) <- "model.info"
  return(model.info)
} )

example(apsrtable)


### Switch back to the default
setMethod("modelInfo", "summary.lm", apsrtable:::modelInfo.summary.lm)
## Not run: 
##D example(apsrtable)
## End(Not run)



cleanEx()
nameEx("notefunctions")
### * notefunctions

flush(stderr()); flush(stdout())

### Name: notefunctions
### Title: Table notes
### Aliases: notefunctions se.note stars.note pval.note

### ** Examples

### Custom note function

signif.pct <- function(env) {
  paste("$^*$ significant at", evalq(lev,envir=env)*100, "percent")
}
### Continue the example from apsrtable
## Not run: 
##D apsrtable(lm.D90, lm.D9, glm.D9, digits=1, align="left",
##D           stars=1, lev=0.05, model.counter=0, order="rl",
##D           notes=list(se.note, signif.pct, 
##D             "Plant weight data from the lm() example" )
##D 	 )
## End(Not run)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
