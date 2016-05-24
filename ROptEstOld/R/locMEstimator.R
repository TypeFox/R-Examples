###############################################################################
## Location M estimator (univariate location)
###############################################################################
setMethod("locMEstimator", signature(x = "numeric", IC = "InfluenceCurve"),
    function(x, IC, eps = .Machine$double.eps^0.5){
        if(numberOfMaps(IC@Curve) > 1)
            stop("number of Maps of 'IC' has to be 1")

        mest <- function(theta, x, IC){
            return(rowSums(evalIC(IC, as.matrix(x-theta))))
        }
        res <- uniroot(f = mest, interval = c(min(x), max(x)), 
             tol = eps, x = x, IC = IC)$root

        return(list(loc = res))
    })
