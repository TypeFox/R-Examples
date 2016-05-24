### R code from vignette source 'RcppClassic-intro.Rnw'

###################################################
### code chunk number 1: RcppClassic-intro.Rnw:36-40
###################################################
require( RcppClassic )
prettyVersion <- packageDescription("RcppClassic")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")
RcppBibfile <- sub("\\.bib$", "", Rcpp:::bib())


###################################################
### code chunk number 3: RcppClassic-intro.Rnw:81-83 (eval = FALSE)
###################################################
## importFrom(Rcpp, evalCpp)
## import(RcppClassic)


