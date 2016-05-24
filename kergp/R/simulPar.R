##' Draw random values for the parameters of a covariance kernel
##' object.


## setMethod("simulPar", 
##           signature = signature(object = "covTS"),
##           definition = function(object, nsim = 1L, seed = NULL){
            
##             ## copied from simulate.lm
##             if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
##               runif(1)                     # initialize the RNG if necessary
##             if(is.null(seed))
##               RNGstate <- get(".Random.seed", envir = .GlobalEnv)
##             else {
##               R.seed <- get(".Random.seed", envir = .GlobalEnv)
##               set.seed(seed)
##               RNGstate <- structure(seed, kind = as.list(RNGkind()))
##               on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
##             }
##             Par <- coef(object)
##             L <- coefLower(object)
##             U <- coefUpper(object)
##             p <- length(Par)
##             Par <- array(NA, dim = c(nsim, p), dimnames = list(NULL, names(Par)))
##             for (k in  1:p) {
##               if ( is.finite(L[k]) ) {
##                 if (is.finite(U[k])) {
##                   Par[ , k] <- runif(n = nsim, min = L[k], max = U[k])
##                 } else {
##                   Par[ , k] <- L[k] + rexp(n = nsim, rate = 1)
##                 } 
##               } else{
##                 if ( is.finite(U[k]) ) {
##                   Par[ , k] <- U[k] - rexp(n = nsim, rate = 1)
##                 } else {
##                   Par[ , k] <- rcauchy(n = nsim)
##                 }   
##               }
##             }
##             Par 
##           }
##           )
