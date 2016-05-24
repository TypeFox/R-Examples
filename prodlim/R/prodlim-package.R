### prodlim-package.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Apr 24 2015 (09:08) 
## Version: 
## last-updated: Dec  1 2015 (11:37) 
##           By: Thomas Alexander Gerds
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Functions for estimating probabilities from right censored data
#'
#' @docType package
#' @name prodlim
#' @useDynLib prodlim
#' @importFrom survival survdiff Surv cluster
#' @importFrom stats quantile 
#' @import lava
#' @importFrom Rcpp sourceCpp
## --> importFrom KernSmooth dpik
#' @importFrom graphics abline axis lines mtext par plot points polygon rect segments strheight strwidth text
#' @importFrom stats .getXlevels delete.response drop.terms formula get_all_vars median model.frame model.matrix model.response na.omit pchisq predict qnorm reformulate terms update update.formula
NULL


#----------------------------------------------------------------------
### prodlim-package.R ends here
