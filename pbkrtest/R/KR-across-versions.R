## #######################################################################################
##
## Functionality required to make pbkrtest work both on CRAN and devel versions of lme4
##
## Banff, August 2013, Søren Højsgaard
##
## #######################################################################################

.get.RT.dim.by.RT <- function(object) {
  ## output: dimension (no of columns) of covariance matrix for random term ii
  .cc <- class(object)
  qq <-
    if (.cc %in% "mer") {
      sapply(object@ST,function(X) nrow(X)) 
    } else {
      sapply(object@cnms, length)  ## FIXME: use getME()
    }
  qq
}

