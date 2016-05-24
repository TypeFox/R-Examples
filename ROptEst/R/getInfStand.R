###############################################################################
## standardizing matrix for asymptotic G-Risk
###############################################################################
setMethod("getInfStand", signature(L2deriv = "UnivariateDistribution",
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, clip, cent, trafo){
        c1 <- cent - clip
        c2 <- cent + clip
        return(trafo/(m2df(L2deriv, c2) - m2df(L2deriv, c1)
                + c1*m1df(L2deriv, c1) - c2*m1df(L2deriv, c2)))
    })
setMethod("getInfStand", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "TotalVarNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, clip, cent, trafo){
        D1 <- sign(as.vector(trafo))*L2deriv
        return(trafo/(m2df(D1, cent+clip) - m2df(D1, cent) + cent*m1df(D1, cent) 
                - (cent+clip)*m1df(D1, cent+clip)))
    })
setMethod("getInfStand", signature(L2deriv = "RealRandVariable",
                                   neighbor = "UncondNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, 
             Distr, A.comp, cent, trafo, w){
        w.fct <- function(x){
            weight(w)(evalRandVar(L2deriv, as.matrix(x)) [,,1]) 
        }
        integrandA <- function(x, L2.i, L2.j, i, j){
            return((L2.i(x) - cent[i])*(L2.j(x) - cent[j])*w.fct(x = x))
        }

        nrvalues <- length(L2deriv)
        erg <- matrix(0, ncol = nrvalues, nrow = nrvalues)
        for(i in 1:nrvalues)
            for(j in i:nrvalues)
                if(A.comp[i,j])
                    erg[i, j] <- E(object = Distr, fun = integrandA, 
                                   L2.i = L2deriv@Map[[i]], 
                                   L2.j = L2deriv@Map[[j]], i = i, j = j, 
                                   useApply = FALSE)

        erg[col(erg) < row(erg)] <- erg[col(erg) > row(erg)]

        return(trafo %*% solve(erg))
    })
###############################################################################
## standardizing constant for one-sided bias
###############################################################################
setMethod("getInfStand", signature(L2deriv = "UnivariateDistribution",
                                   neighbor = "ContNeighborhood",
                                   biastype = "onesidedBias"),
    function(L2deriv, neighbor, biastype, clip, cent, trafo){
        c1 <- if (sign(biastype)<0) cent - clip else -Inf
        c2 <- if (sign(biastype)>0) cent + clip else Inf
        m1 <- if (sign(biastype)<0) m2df(L2deriv, c1) else 0
        m2 <- if (sign(biastype)>0) m2df(L2deriv, c2) else E(L2deriv, function(x)x^2)
        c10 <- if (sign(biastype)<0) c1*m1df(L2deriv, c1) else 0
        c20 <- if (sign(biastype)>0) c2*m1df(L2deriv, c2) else 0
        return(trafo/(m2 - m1 + c10 - c20))
    })

###############################################################################
## standardizing constant for asymmetric bias
###############################################################################
setMethod("getInfStand", signature(L2deriv = "UnivariateDistribution",
                                   neighbor = "ContNeighborhood",
                                   biastype = "asymmetricBias"),
    function(L2deriv, neighbor, biastype, clip, cent, trafo){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]
        c1 <- cent - clip/nu1
        c2 <- cent + clip/nu2
        return(trafo/(m2df(L2deriv, c2) - m2df(L2deriv, c1)
                + c1*m1df(L2deriv, c1) - c2*m1df(L2deriv, c2)))
    })
