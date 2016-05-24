###############################################################################
## centering constant for asymptotic MSE and asymptotic Hampel risk 
###############################################################################
setMethod("getInfCent", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "ContNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, 
             clip, cent, tol.z, symm, trafo){
        if(symm) return(0)

        z.fct <- function(z, c0, D1){
            return(c0 + (z-c0)*p(D1)(z-c0) - (z+c0)*p(D1)(z+c0) + m1df(D1, z+c0) - m1df(D1, z-c0))
        }
        lower <- getLow(L2deriv)
        upper <- getUp(L2deriv)

        return(uniroot(z.fct, lower = lower, upper = upper, tol = tol.z, 
                    c0=clip, D1=L2deriv)$root)
    })
setMethod("getInfCent", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "TotalVarNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, 
             clip, cent, tol.z, symm, trafo){
        if(symm) return(-clip/2)

        D1 <- sign(as.vector(trafo))*L2deriv
        g.fct <- function(g, c0, D1){
            return(g*p(D1)(g) + (g+c0)*(p(D1)(g+c0, lower.tail = FALSE)) - m1df(D1, g) + m1df(D1, g+c0))
        }
        lower <- -clip
        upper <- 0
        return(uniroot(g.fct, lower = lower, upper = upper, tol = tol.z,
                    c0 = clip, D1 = D1)$root)
    })


setMethod("getInfCent", signature(L2deriv = "RealRandVariable",
                                  neighbor = "TotalVarNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, Distr, z.comp, w,
             tol.z = .Machine$double.eps^.5){
        stand <- stand(w)
        clip <- clip(w)
        b <- clip[2]-clip[1]
        ### if(symm) return(b/2)

        g.fct <- function(g, c0){
            fct <- function(x){
                  Lx <- evalRandVar(L2deriv, as.matrix(x)) [,,1]
                  Y <- as.numeric(stand%*%Lx)
                  pmin(pmax(g,Y),g+c0)
                  }
            return(E(object = Distr, fun = fct, useApply = FALSE))
        }
        lower <- -b
        upper <- 0

        return(uniroot(g.fct, lower = lower, upper = upper, tol = tol.z,
                    c0 = b)$root)
    })

setMethod("getInfCent", signature(L2deriv = "RealRandVariable",
                                  neighbor = "ContNeighborhood",
                                  biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, Distr, z.comp, w,
             tol.z = .Machine$double.eps^.5){
        integrand1 <- function(x){
            weight(w)(evalRandVar(L2deriv, as.matrix(x)) [,,1]) 
        }
        integrand2 <- function(x, L2.i){
            return(L2.i(x)*integrand1(x))
        }

        res1 <- E(object = Distr, fun = integrand1, useApply = FALSE)
        nrvalues <- length(L2deriv)
        res2 <- numeric(nrvalues)
        for(i in 1:nrvalues){
            if(z.comp[i]){
                 res2[i] <- E(object = Distr, fun = integrand2, 
                              L2.i = L2deriv@Map[[i]], useApply = FALSE)
            }else{            
                res2[i] <- 0
            }
        }

        return(res2/res1)
    })
###############################################################################
## centering constant for asymptotic one-sided MSE and asymptotic one-sided Hampel risk
###############################################################################
setMethod("getInfCent", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "ContNeighborhood",
                                  biastype = "onesidedBias"),
    function(L2deriv, neighbor, biastype, clip, cent, tol.z, symm, trafo){
        if (sign(biastype)> 0){
        z.fct <- function(z, c0, D1){
            return(c0 - (z+c0)*p(D1)(z+c0) + m1df(D1, z+c0))
        }
        lower <- getLow(L2deriv)
        upper <- 0
        }else{
        z.fct <- function(z, c0, D1){
            return(- z + (z-c0)*p(D1)(z-c0) - m1df(D1, z-c0))
        }
        lower <- 0
        upper <- getUp(L2deriv)
        }
        return(uniroot(z.fct, lower = lower, upper = upper, tol = tol.z,
                    c0=clip, D1=L2deriv)$root)
    })

setMethod("getInfCent", signature(L2deriv = "UnivariateDistribution",
                                  neighbor = "ContNeighborhood",
                                  biastype = "asymmetricBias"),
    function(L2deriv, neighbor, biastype, clip, cent, tol.z, symm, trafo){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]

        z.fct <- function(z, c0, D1){
            return(c0/nu2 + (z-c0/nu1)*p(D1)(z-c0/nu1) -
                   (z+c0/nu2)*p(D1)(z+c0/nu2) + m1df(D1, z+c0/nu2) -
                   m1df(D1, z-c0/nu1))
        }
        lower <- getLow(L2deriv)
        upper <- getUp(L2deriv)

        return(uniroot(z.fct, lower = lower, upper = upper, tol = tol.z,
                    c0=clip, D1=L2deriv)$root)
    })
