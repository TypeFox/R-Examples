###############################################################################
## asyVar 
###############################################################################
setMethod("getInfV", signature(L2deriv = "UnivariateDistribution",
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, clip, cent, stand){
        c1 <- cent - clip
        c2 <- cent + clip
        return(stand^2*(m2df(L2deriv, c2) - m2df(L2deriv, c1)
                + 2 * cent *(m1df(L2deriv, c1) - m1df(L2deriv, c2))
                + cent^2 * (p(L2deriv)(c2) -p(L2deriv)(c1))
                + clip^2 * (p(L2deriv)(c2, lower.tail = FALSE) +p(L2deriv)(c1))
                ))
    })


setMethod("getInfV", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "TotalVarNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, clip, cent, stand){
        c1 <- cent
        c2 <- cent+clip
        return(stand^2*(m2df(L2deriv, c2) - m2df(L2deriv, c1)
                + c2^2 * (p(L2deriv)(c2, lower.tail = FALSE))
                + c1^2* p(L2deriv)(c1)
                ))
    })

setMethod("getInfV", signature(L2deriv = "RealRandVariable",
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, Distr, V.comp, 
             cent, stand, w){
        w.fct <- function(x){
            (weight(w)(evalRandVar(L2deriv, as.matrix(x)) [,,1]))^2 
        }
        
        cent0 <- solve(stand, cent)


        integrandV <- function(x, L2.i, L2.j, i, j){
            return((L2.i(x) - cent0[i])*(L2.j(x) - cent0[j])*w.fct(x = x))
        }

        nrvalues <- length(L2deriv)
        erg <- matrix(0, ncol = nrvalues, nrow = nrvalues)
        for(i in 1:nrvalues)
            for(j in i:nrvalues)
                if(V.comp[i,j])
                    erg[i, j] <- E(object = Distr, fun = integrandV, 
                                   L2.i = L2deriv@Map[[i]], 
                                   L2.j = L2deriv@Map[[j]], i = i, j = j, 
                                   useApply = FALSE)

        erg[col(erg) < row(erg)] <- erg[col(erg) > row(erg)]

        return(as.matrix(stand %*% erg %*% t(stand)))
    })
setMethod("getInfV", signature(L2deriv = "RealRandVariable",
                                   neighbor = "TotalVarNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, Distr, V.comp,
             cent, stand, w){
        w.fct <- function(x){
            (weight(w)(evalRandVar(L2deriv, as.matrix(x)) [,,1]))^2
        }


        integrandV <- function(x){
            L2 <- evalRandVar(L2deriv, as.matrix(x)) [,,1]
            Y <- stand %*% L2
            return(Y^2 * w.fct(x = x))
        }

        return(matrix(E(object = Distr, fun = integrandV, useApply = FALSE),
                        ncol = 1, nrow = 1))

    })
###############################################################################
## standardizing constant for one-sided bias
###############################################################################
setMethod("getInfV", signature(L2deriv = "UnivariateDistribution",
                                   neighbor = "ContNeighborhood",
                                   biastype = "onesidedBias"),
    function(L2deriv, neighbor, biastype, clip, cent, stand){
        c1 <- if (sign(biastype)<0) cent - clip else -Inf
        c2 <- if (sign(biastype)>0) cent + clip else Inf
        V1 <- if (sign(biastype)<0) m2df(L2deriv, c1) else 0
        V2 <- if (sign(biastype)>0) m2df(L2deriv, c2) else E(L2deriv, function(x)x^2)
        E1 <- if (sign(biastype)<0) m2df(L2deriv, c1) else 0
        E2 <- if (sign(biastype)>0) m2df(L2deriv, c2) else 0
        F1 <- if (sign(biastype)<0) p(L2deriv)(c1) else 0
        F2 <- if (sign(biastype)>0) p(L2deriv)(c2, lower.tail = FALSE) else 0
        c10 <- if (sign(biastype)<0) c1*m1df(L2deriv, c1) else 0
        c20 <- if (sign(biastype)>0) c2*m1df(L2deriv, c2) else 0
        c1 <- cent - clip
        c2 <- cent + clip
        return(stand^2*(V2 - V1 + 2 * cent * (E1 - E2) 
                + cent^2 * (1 - F2 - F1)
                + clip^2 * (F2 + F1)
                ))
    })


###############################################################################
## standardizing constant for asymmetric bias
###############################################################################
setMethod("getInfV", signature(L2deriv = "UnivariateDistribution",
                                   neighbor = "ContNeighborhood",
                                   biastype = "asymmetricBias"),
    function(L2deriv, neighbor, biastype, clip, cent, stand){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]
        c1 <- cent - clip/nu1
        c2 <- cent + clip/nu2
        V0 <- m2df(L2deriv, c2) - m2df(L2deriv, c1)
        V1 <- m1df(L2deriv, c2) - m1df(L2deriv, c1)
        V2 <- p(L2deriv)(c2) -p(L2deriv)(c1)
        V3 <- (p(L2deriv)(c2, lower.tail=FALSE))/nu2^2 +p(L2deriv)(c1)/nu1^2
        V <- stand^2*( V0 - 2 * cent * V1 + cent^2 * V2 + clip^2 * V3)
        return(V)
    })
