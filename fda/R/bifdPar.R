bifdPar = function(bifdobj, Lfdobjs=int2Lfd(2), Lfdobjt=int2Lfd(2), 
                   lambdas=0, lambdat=0, estimate=TRUE) {
# Sets up a bivariate functional parameter object
#  Arguments:
#  BIFDOBJ  ... A bivariate functional data object.  The basis for this object 
#               is used to define the bivariate functional parameter.
#               When an initial value is required for iterative 
#               estimation of a bivariate functional parameter, the coefficients
#               will give the initial values for the iteration.
#  LFDOBJS  ... A linear differential operator value or a derivative
#               value for penalizing the roughness of the object
#               with respect to the first argument s.
#               By default, this is 2.
#  LFDOBJT  ... A linear differential operator value or a derivative
#               value for penalizing the roughness of the object
#               with respect to the second argument t.
#               By default, this is 2.
#  LAMBDAS  ... The penalty parameter controlling the smoothness of
#               the estimated parameter with respect to the first argument s.  
#               By default this is 0.
#  LAMBDAT  ... The penalty parameter controlling the smoothness of
#               the estimated parameter with respect to the second argument t.  
#               By default this is 0.
#  ESTIMATE ... If nonzero, the parameter is estimated; if zero, the
#                parameter is held fixed at this value.
#                By default, this is 1.

#  last modified 28 October 2009

#  check BIFDOBJ

if (!inherits(bifdobj, "bifd")) {
    stop("BIFDOBJ is not a bivariate functional data object.")
}

#  check the linear differential operators

Lfdobjs = int2Lfd(Lfdobjs)
Lfdobjt = int2Lfd(Lfdobjt)

if (!is.Lfd(Lfdobjs)) {
    stop("LFDOBJS is not a linear differential operator object.")
}
if (!is.Lfd(Lfdobjt)) {
    stop("LFDOBJT is not a linear differential operator object.")
}

#  check the roughness penalty parameters

if (!is.numeric(lambdas) ) {
    stop("LAMBDAS is not numeric.")
}
if (lambdas < 0) {
    warning("LAMBDAS is negative, and is set to zero.")
    lambdas = 0
}
if (!is.numeric(lambdat)) {
    stop("LAMBDAT is not numeric.")
}
if (lambdat < 0) {
    warning("LAMBDAT is negative, and is set to zero.")
    lambdat = 0
}

if (!is.logical(estimate)) {
    stop("ESTIMATE is not logical.")
}

#  set up the bifdPar object

bifdParobj <- list(bifd=bifdobj, estimate=estimate, 
                   lambdas=lambdas, lambdat=lambdat, Lfds=Lfdobjs, Lfdt=Lfdobjt, 
                   estimate=estimate)

oldClass(bifdParobj) <- "bifdPar"

return(bifdParobj)

}
