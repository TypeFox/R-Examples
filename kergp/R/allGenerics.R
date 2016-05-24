##***********************************************************************
## Some new generics seem to be required to reach portability.
## 
##***********************************************************************


##=======================================================================
## Map the parameter vector with the kernel parameters.
##===========================================================A============
setGeneric("parMap",
           function(object, ...) standardGeneric("parMap")
           )
setGeneric("npar",
           function(object, ...) standardGeneric("npar")
           )

##=======================================================================
## Check that the design is compatible with the object based on names,
## dim, etc.
##=======================================================================
setGeneric("checkX",
           function(object, X, ...) standardGeneric("checkX")
           )

##=======================================================================
## Extract the official name of (USUALLY ONE-DIMENSIONAL) kernel. 
##=======================================================================
if (!isGeneric("kernelName")) {
  setGeneric("kernelName",
             function(object, ...) standardGeneric("kernelName")
             )
}

##=======================================================================
## Extract the names of the inputs. The are usesd for 'official'
## validations in kergp.
##=======================================================================
if (!isGeneric("inputNames")) {
  setGeneric("inputNames",
             function(object, ...) standardGeneric("inputNames")
             )
}

##=======================================================================
## replacement method for coef
##=======================================================================
if (!isGeneric("coef<-")) {
  setGeneric("coef<-",
           function(object, ..., value) standardGeneric("coef<-")
           )
}
##=======================================================================
## Extract or set bounds on parameters
##=======================================================================
setGeneric("coefLower",
           function(object, ...) standardGeneric("coefLower")
           )
setGeneric("coefLower<-",
           function(object, ..., value) standardGeneric("coefLower<-")
           )
setGeneric("coefUpper",
           function(object, ...) standardGeneric("coefUpper")
           )
setGeneric("coefUpper<-",
           function(object, ..., value) standardGeneric("coefUpper<-")
           )

##=======================================================================
## Compute bounds for parameters from a given design (identifiability)
##
## XXX the name is not very appealing... this could be simply coefLower
## with amultiple dispatch on 'object' and 'X'. Then change the signature
## of the generic 'coefLower' and 'coefUpper' to inclue an 'X' formal. 
##
##=======================================================================
setGeneric("compCoefLower",
           function(object, X, ...) standardGeneric("compCoefLower")
           )
setGeneric("compCoefUpper",
           function(object, X, ...) standardGeneric("compCoefUpper")
           )

##=======================================================================
## The same as covMatrix, but with more arguements
##=======================================================================
setGeneric("covMat",
           function(object, X, Xnew = NULL, ...) standardGeneric("covMat")
           )

##=======================================================================
## draw parameters at random from a covariance structure
##=======================================================================
setGeneric("simulPar",
           function(object, nsim = 1L, seed = NULL, ...) standardGeneric("simulPar")
           )

##=======================================================================
## gls fit from a covariance object. Works for an instance of covariance
## kernel, but could be made to work on a matrix.
##=======================================================================
setGeneric("gls",
           function(object, ...) standardGeneric("gls")
           )

##==================================
## mle fit from a covariance object. 
##==================================
if (!isGeneric("mle")) {
  setGeneric("mle",
             function(object, ...) standardGeneric("mle")
  )
}

##==================================
## scores 
##==================================
if (!isGeneric("scores")) {
  setGeneric("scores",
             function(object, ...) standardGeneric("scores")
             )
}




