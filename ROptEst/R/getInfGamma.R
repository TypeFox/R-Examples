###############################################################################
## gamma in case of a convex asymptotic risk
###############################################################################
setMethod("getInfGamma", signature(L2deriv = "UnivariateDistribution",
                                   risk = "asGRisk", 
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, risk, neighbor, biastype, cent, clip){
        c1 <- cent - clip
        c2 <- cent + clip
        return(m1df(L2deriv, c2) + m1df(L2deriv, c1) 
                    - c1*p(L2deriv)(c1) + c2*p(L2deriv)(c2, lower.tail = FALSE))
    })
###############################################################################
## r^2 b = E(c - A Lambda)_+ Probleme mit Startwerten!!!
## daher: r^2 b = E(A Lambda - (c+b))_+ 
###############################################################################
setMethod("getInfGamma", signature(L2deriv = "UnivariateDistribution",
                                   risk = "asGRisk", 
                                   neighbor = "TotalVarNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, risk, neighbor, biastype, cent, clip){
        return(m1df(L2deriv, cent+clip) + (cent+clip)*p(L2deriv)(cent+clip,
               lower.tail = FALSE))
    })

setMethod("getInfGamma", signature(L2deriv = "RealRandVariable",
                                   risk = "asMSE", 
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, risk, neighbor, biastype, Distr, 
             stand, cent, clip, power = 1L){
        integrandG <- function(x, L2, stand, cent, clip){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- stand %*% X
            res <- norm(risk)(Y) - clip

            return((res > 0)*res^power)
        }

        return(-E(object = Distr, fun = integrandG, L2 = L2deriv, 
                  stand = stand, cent = cent, clip = clip, useApply = FALSE))
    })

setMethod("getInfGamma", signature(L2deriv = "RealRandVariable",
                                   risk = "asMSE",
                                   neighbor = "TotalVarNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, risk, neighbor, biastype, Distr,
             stand, cent, clip, power = 1L){
        integrandG <- function(x, L2, stand, cent, clip){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- stand %*% X
            res <- Y - clip

            return((res > 0)*res^power)
        }

        return(-E(object = Distr, fun = integrandG, L2 = L2deriv,
                  stand = stand, cent = cent, clip = clip, useApply = FALSE))
    })
###############################################################################
## gamma in case of asymptotic under-/overshoot risk
###############################################################################
setMethod("getInfGamma", signature(L2deriv = "UnivariateDistribution",
                                   risk = "asUnOvShoot", 
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, risk, neighbor, biastype, cent, clip){
        return(2*(m1df(L2deriv, cent+clip) + (cent+clip)*(1-p(L2deriv)(cent+clip))))
    })

###############################################################################
## gamma in case of asymptotic one-sided convex asymptotic risk
###############################################################################
setMethod("getInfGamma", signature(L2deriv = "UnivariateDistribution",
                                   risk = "asMSE",
                                   neighbor = "ContNeighborhood",
                                   biastype = "onesidedBias"),
    function(L2deriv, risk, neighbor, biastype, cent, clip){
        c1 <- cent - clip 
        c2 <- cent + clip 
        if (sign(biastype)<0) 
           return (m1df(L2deriv, c1) -c1*p(L2deriv)(c1))
        else 
           return (m1df(L2deriv, c2) +c2*(1-p(L2deriv)(c2)))
    })

###############################################################################
## gamma in case of a asymmetric asymptotic risk
###############################################################################
setMethod("getInfGamma", signature(L2deriv = "UnivariateDistribution",
                                   risk = "asMSE",
                                   neighbor = "ContNeighborhood",
                                   biastype = "asymmetricBias"),
    function(L2deriv, risk, neighbor, biastype, cent, clip){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]

        c1 <- cent - clip/nu1
        c2 <- cent + clip/nu2
        return(m1df(L2deriv, c2)/nu2 + m1df(L2deriv, c1)/nu1
                    - c1*p(L2deriv)(c1)/nu1 + c2*(1-p(L2deriv)(c2))/nu2)
    })

