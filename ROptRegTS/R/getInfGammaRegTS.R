###############################################################################
## gamma in case of a convex asymptotic risk and 
## (unconditional) contamination neighborhoods (* = c)
###############################################################################
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution",
                                        risk = "asMSE", 
                                        neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip){
        print("falsch")
        Gfct <- function(x, cent, clip, D1){
            if(x == 0) return(max(abs(cent) - clip, 0))
            c1 <- cent/x - clip/abs(x)
            c2 <- cent/x + clip/abs(x)
                        
            return(abs(x)*(m1df(D1, c2) + m1df(D1, c1) - c1*p(D1)(c1) + c2*(1-p(D1)(c2))))
        }

        return(E(Regressor, Gfct, cent = cent, clip = clip, D1 = ErrorL2deriv))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "MultivariateDistribution",
                                        risk = "asMSE", 
                                        neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip){
        if(z.comp){
            Gfct <- function(x, stand, cent, clip, D1){
                Gfctu <- function(u, xx, stand, cent, clip){
                    v <- as.vector(stand %*% (xx*u - cent))
                    res <- as.vector(sqrt(v %*% v)) - clip
                    return((res > 0)*res)
                }
                E(D1, Gfctu, xx = x, stand = stand, cent = cent, clip = clip)
            }
            return(-E(Regressor, Gfct, stand = stand, cent = cent, clip = clip, D1 = ErrorL2deriv))
        }else{
            Gfct <- function(x, stand, clip, D1){
                v <- t(x) %*% stand
                v <- as.vector(sqrt(v %*% t(v)))
                c0 <- clip/v

                return(v*(m1df(D1, c0) + m1df(D1, -c0) + 2*c0*p(D1)(-c0)))
            }
            return(E(Regressor, Gfct, stand = stand, clip = clip, D1 = ErrorL2deriv))
        }
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution",
                                        risk = "asMSE", 
                                        neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip){
        Gfct <- function(x, cent, clip, D1){
            if(x == 0) return(max(abs(cent(x)) - clip, 0))
            c1 <- cent(x) - clip/abs(x)
            c2 <- cent(x) + clip/abs(x)

            return(abs(x)*(m1df(D1, c2) + m1df(D1, c1) - c1*p(D1)(c1) + c2*(1-p(D1)(c2))))
        }
        return(E(Regressor, Gfct, cent = cent, clip = clip, D1 = ErrorL2deriv))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "MultivariateDistribution",
                                        risk = "asMSE", 
                                        neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip){
        if(z.comp){
            Gfct <- function(x, stand, cent, clip, D1){
                Gfctu <- function(u, xx, stand, cent, clip){
                    v <- as.vector(stand %*% xx*(u - cent(xx)))
                    res <- as.vector(sqrt(v %*% v)) - clip
                    return((res > 0)*res)
                }
                E(D1, Gfctu, xx = x, stand = stand, cent = cent, clip = clip)
            }
            return(-E(Regressor, Gfct, stand = stand, cent = cent, clip = clip, D1 = ErrorL2deriv))
        }else{
            Gfct <- function(x, stand, clip, D1){
                v <- t(x) %*% stand
                v <- as.vector(sqrt(v %*% t(v)))
                c0 <- clip/v

                return(v*(m1df(D1, c0) + m1df(D1, -c0) + 2*c0*p(D1)(-c0)))
            }
            return(E(Regressor, Gfct, stand = stand, clip = clip, D1 = ErrorL2deriv))
        }
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "Distribution",
                                        risk = "asMSE", 
                                        neighbor = "Av2CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip){
        c1 <- cent - clip
        c2 <- cent + clip
        return(m1df(ErrorL2deriv, c2) + m1df(ErrorL2deriv, c1) 
                    - c1*p(ErrorL2deriv)(c1) + c2*(1-p(ErrorL2deriv)(c2)))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution",
                                        risk = "asMSE", 
                                        neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip){
        Gfct <- function(x, cent, clip, D1){
            if(x == 0) return(0)

            c0 <- (cent(x) + clip)/abs(x)

            return((m1df(D1, c0) + c0*(1-p(D1)(c0)))*abs(x))
        }
        return(E(Regressor, Gfct, cent = cent, clip = clip, D1 = ErrorL2deriv))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "MultivariateDistribution",
                                        risk = "asMSE", 
                                        neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent, clip){
        Gfct <- function(x, cent, clip, stand, D1){
            if(all(x == numeric(length(x)))) return(0)

            v <- as.vector(sqrt(sum((stand %*% x)^2)))
            c0 <- (cent(x) + clip)/v
            
            return((m1df(D1, c0) + c0*(1-p(D1)(c0)))*v)
        }
        return(E(Regressor, Gfct, cent = cent, clip = clip, stand = stand, D1 = ErrorL2deriv))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "RealRandVariable",
                                        Regressor = "Distribution",
                                        risk = "asMSE", 
                                        neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorDistr, stand, cent, clip){
        integrandGu <- function(u, xx, ErrorL2deriv, A, z, b, k0){ 
            L1 <- as.vector(evalRandVar(ErrorL2deriv, u))
            L <- c(xx*L1[1], L1[2:k0])
            Y <- A %*% (L - z)
            res <- as.vector(sqrt(sum(Y^2))) - b

            return((res > 0)*res)
        }
        integrandGx <- function(x, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandGu, xx = x, 
              ErrorL2deriv = ErrorL2deriv, A = stand, z = cent, b = clip, k0 = k0)
        }
        k <- length(ErrorL2deriv)
    
        return(-E(object = Regressor, fun = integrandGx, ErrorL2deriv = ErrorL2deriv, 
                  ErrorDistr = ErrorDistr, stand = stand, cent = cent, clip = clip, k0 = k))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "RealRandVariable",
                                        Regressor = "Distribution",
                                        risk = "asMSE", 
                                        neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorDistr, stand, cent, clip){
        k <- length(ErrorL2deriv)
        integrandGu <- function(u, xx, ErrorL2deriv, A, z, b, k0){ 
            L1 <- as.vector(evalRandVar(ErrorL2deriv, u))
            z1 <- z(xx)
            L <- c(xx*L1[1]-z1[1], L1[2:k0]-z1[2:k0])
            Y <- as.vector(A %*% L)
            res <- as.vector(sqrt(sum(Y^2))) - b

            return((res > 0)*res)
        }
        integrandGx <- function(x, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandGu, xx = x, 
              ErrorL2deriv = ErrorL2deriv, A = stand, z = cent, b = clip, k0 = k0)
        }

        return(-E(object = Regressor, fun = integrandGx, ErrorL2deriv = ErrorL2deriv, 
                  ErrorDistr = ErrorDistr, stand = stand, cent = cent, clip = clip, k0 = k))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution",
                                        risk = "asUnOvShoot", 
                                        neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, cent, clip){
        Gfct <- function(x, cent, clip, D1){
            if(x == 0) return(max(-(cent + clip), 0))
            c1 <- (cent + clip)/x
            return(max(x,0)*(m1df(D1, c1) + c1*(1-p(D1)(c1)))
                   - min(x, 0)*(m1df(D1, c1) - c1*p(D1)(c1)))
        }
        return(2*E(Regressor, Gfct, cent = cent, clip = clip, D1 = ErrorL2deriv))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution",
                                        risk = "asUnOvShoot", 
                                        neighbor = "TotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, cent, clip){
        Gfct <- function(x, cent, clip, D1){
            if(x == 0) return(0)
            c1 <- (cent + clip)/x
            return(max(x,0)*(m1df(D1, c1) + c1*(1-p(D1)(c1)))
                   - min(x, 0)*(m1df(D1, c1) - c1*p(D1)(c1)))
        }
        return(E(Regressor, Gfct, cent = cent, clip = clip, D1 = ErrorL2deriv))
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "numeric",
                                        risk = "asUnOvShoot", 
                                        neighbor = "CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, clip){
        bx <- clip/Regressor
        if(Regressor > 0){
            return(2*(Regressor*m1df(ErrorL2deriv, bx) + clip*(1 - p(ErrorL2deriv)(bx))))
        }else{
            return(-2*(Regressor*m1df(ErrorL2deriv, bx) - clip*p(ErrorL2deriv)(bx)))
        }
    })
setMethod("getInfGammaRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "numeric",
                                        risk = "asUnOvShoot", 
                                        neighbor = "CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, clip){
        bx <- clip/Regressor
        if(Regressor > 0){
            return(Regressor*m1df(ErrorL2deriv, bx) + clip*(1 - p(ErrorL2deriv)(bx)))
        }else{
            return(-(Regressor*m1df(ErrorL2deriv, bx) - clip*p(ErrorL2deriv)(bx)))
        }
    })
