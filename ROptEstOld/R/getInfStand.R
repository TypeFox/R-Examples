###############################################################################
## standardizing matrix for asymptotic G-Risk
###############################################################################
setMethod("getInfStand", signature(L2deriv = "UnivariateDistribution",
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, neighbor, clip, cent, trafo){
        c1 <- cent - clip
        c2 <- cent + clip
        return(trafo/(m2df(L2deriv, c2) - m2df(L2deriv, c1)
                + c1*m1df(L2deriv, c1) - c2*m1df(L2deriv, c2)))
    })
setMethod("getInfStand", signature(L2deriv = "UnivariateDistribution",
                                neighbor = "TotalVarNeighborhood"),
    function(L2deriv, neighbor, clip, cent, trafo){
        D1 <- sign(as.vector(trafo))*L2deriv
        return(trafo/(m2df(D1, cent+clip) - m2df(D1, cent) + cent*m1df(D1, cent) 
                - (cent+clip)*m1df(D1, cent+clip)))
    })
setMethod("getInfStand", signature(L2deriv = "RealRandVariable",
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, neighbor, Distr, A.comp, stand, clip, cent, trafo){
        w.fct <- function(x, L2, stand, cent, clip){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- apply(X, 2, "%*%", t(stand)) 
            h.vct <- sqrt(colSums(Y^2))
            ind2 <- (h.vct < clip/2)
            h.vct <- ind2*clip/2 + (1-ind2)*h.vct
            ind1 <- (h.vct < clip)

            return(ind1 + (1-ind1)*clip/h.vct)
        }
        integrandA <- function(x, L2.i, L2.j, i, j, L2, stand, cent, clip){
            return((L2.i(x) - cent[i])*(L2.j(x) - cent[j])*w.fct(x = x, L2 = L2, 
                                                stand = stand, cent = cent, clip = clip))
        }

        nrvalues <- length(L2deriv)
        erg <- matrix(0, ncol = nrvalues, nrow = nrvalues)
        for(i in 1:nrvalues)
            for(j in i:nrvalues)
                if(A.comp[i,j])
                    erg[i, j] <- E(object = Distr, fun = integrandA, L2.i = L2deriv@Map[[i]], 
                                   L2.j = L2deriv@Map[[j]], i = i, j = j, L2 = L2deriv, 
                                   stand = stand, cent = cent, clip = clip, 
                                   useApply = FALSE)

        erg[col(erg) < row(erg)] <- erg[col(erg) > row(erg)]

        return(trafo %*% solve(erg))
    })
