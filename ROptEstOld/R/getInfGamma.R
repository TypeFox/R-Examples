###############################################################################
## gamma in case of a convex asymptotic risk
###############################################################################
setMethod("getInfGamma", signature(L2deriv = "UnivariateDistribution",
                                   risk = "asMSE", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, cent, clip){
        c1 <- cent - clip
        c2 <- cent + clip
        return(m1df(L2deriv, c2) + m1df(L2deriv, c1) 
                    - c1*p(L2deriv)(c1) + c2*(1-p(L2deriv)(c2)))
    })
###############################################################################
## r^2 b = E(c - A Lambda)_+ Probleme mit Startwerten!!!
## daher: r^2 b = E(A Lambda - (c+b))_+ 
###############################################################################
setMethod("getInfGamma", signature(L2deriv = "UnivariateDistribution",
                                   risk = "asGRisk", 
                                   neighbor = "TotalVarNeighborhood"),
    function(L2deriv, risk, neighbor, cent, clip){
        return(m1df(L2deriv, cent+clip) + (cent+clip)*(1-p(L2deriv)(cent+clip)))
    })
setMethod("getInfGamma", signature(L2deriv = "RealRandVariable",
                                   risk = "asMSE", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, stand, cent, clip){
        integrandG <- function(x, L2, stand, cent, clip){ 
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- apply(X, 2, "%*%", t(stand)) 
            res <- sqrt(colSums(Y^2)) - clip

            return((res > 0)*res)
        }

        return(-E(object = Distr, fun = integrandG, L2 = L2deriv, 
                  stand = stand, cent = cent, clip = clip, useApply = FALSE))
    })

###############################################################################
## gamma in case of asymptotic under-/overshoot risk
###############################################################################
setMethod("getInfGamma", signature(L2deriv = "UnivariateDistribution",
                                   risk = "asUnOvShoot", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, cent, clip){
        return(2*(m1df(L2deriv, cent+clip) + (cent+clip)*(1-p(L2deriv)(cent+clip))))
    })
