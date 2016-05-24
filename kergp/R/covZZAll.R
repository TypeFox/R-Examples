## ------------
## Union Class
## ------------

setClassUnion("covAll", c("covTS", "covMan"))


##***********************************************************************
## 'inputNames'
##
## XXX   This method should be removed for classes inheriting from
##       "covKernel", excpet if inputnames are stored under another
##       name.
## 
##***********************************************************************
setMethod("inputNames",
          signature = signature(object = "covAll"),
          definition = function(object, ...){
            object@inputNames
          })




##***********************************************************************
## Check that the design matrix X is compatible with the covariance
## object. 
##***********************************************************************
setMethod("checkX",
          signature = signature(object = "covAll", X = "matrix"),
          definition = function(object, X, strict = FALSE, ...){
            iN <- inputNames(object)
            if (strict) {
              if (length(iN) != ncol(X) || !all(iN == colnames(X)) )
                stop("colnames(X) and inpuNames(object) must be identical")
            }
            if ( !all(inputNames(object) %in% colnames(X)) )
              stop("all the elements of inputNames(object) must be in colnames(X)")
            X[ , iN, drop = FALSE] 
          })


##' Draw random values for the parameters of a covariance kernel
##' object.
##'
##' Draw random values for the parameters of a covariance kernel
##' object using the informations \code{coefLower} and
##' \code{coefUpper}.
##' @title Draw random values for the parameters of a covariance kernel.
##' @param object A covariance kernel.
##' @param nsim Number of drawings.
##' @param seed Seed for the random generator.
##' @param ... Other arguments for methods.
##' @return A matrix with \code{nsim} rows and \code{npar(object)} columns.
##' 
##' @author ODY

setMethod("simulPar", 
          signature = signature(object = "covAll"),
          definition = function(object, nsim = 1L, seed = NULL){
            
            ## copied from simulate.lm
            if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
              runif(1)                     # initialize the RNG if necessary
            if(is.null(seed))
              RNGstate <- get(".Random.seed", envir = .GlobalEnv)
            else {
              R.seed <- get(".Random.seed", envir = .GlobalEnv)
              set.seed(seed)
              RNGstate <- structure(seed, kind = as.list(RNGkind()))
              on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
            }
            Par <- coef(object)
            L <- coefLower(object)
            U <- coefUpper(object)
            p <- length(Par)
            Par <- array(NA, dim = c(nsim, p), dimnames = list(NULL, names(Par)))
            for (k in  1:p) {
              if ( is.finite(L[k]) ) {
                if (is.finite(U[k])) {
                  Par[ , k] <- runif(n = nsim, min = L[k], max = U[k])
                } else {
                  Par[ , k] <- L[k] + rexp(n = nsim, rate = 1)
                } 
              } else{
                if ( is.finite(U[k]) ) {
                  Par[ , k] <- U[k] - rexp(n = nsim, rate = 1)
                } else {
                  Par[ , k] <- rcauchy(n = nsim)
                }   
              }
            }
            Par 
          }
          )







