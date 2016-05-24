###############################################################################
## standardizing matrix for asymptotic MSE (* = c)
###############################################################################
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution", 
                                        neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, z.comp, clip, cent, stand, trafo){
        Afct <- function(x, cent, clip, D1, trafo){
            if(x == 0) return(0)
            c1 <- cent/x - clip/abs(x)
            c2 <- cent/x + clip/abs(x)

            return(x^2*(m2df(D1, c2) - m2df(D1, c1) + c1*m1df(D1, c1) - c2*m1df(D1, c2)))
        }
        return(as.matrix(trafo/E(Regressor, Afct, cent = cent, clip = clip, D1 = ErrorL2deriv,
                 trafo = trafo)))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution", 
                                        neighbor = "TotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent){
        Afct <- function(x, cent, clip, D1){
            if(x == 0) return(0)
            c1 <- cent/x
            c2 <- (cent + clip)/x
            res1 <- max(x, 0)^2*(c1*m1df(D1, c1) + m2df(D1, c2) - m2df(D1, c1) - c2*m1df(D1, c2))
            res2 <- min(x, 0)^2*(c2*m1df(D1, c2) + m2df(D1, c1) - m2df(D1, c2) - c1*m1df(D1, c1))

            return(res1 + res2)
        }

        return(1/E(Regressor, Afct, cent = cent, clip = clip, D1 = ErrorL2deriv))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution", 
                                        neighbor = "CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent){
        Afct <- function(x, cent, clip, D1){
            if(x == 0) return(0)
            c1 <- cent(x)/x
            c2 <- clip(x)/x
            res1 <- max(x, 0)^2*(c1*m1df(D1, c1) + m2df(D1, c2) - m2df(D1, c1) - c2*m1df(D1, c2))
            res2 <- min(x, 0)^2*(c2*m1df(D1, c2) + m2df(D1, c1) - m2df(D1, c2) - c1*m1df(D1, c1))

            return(res1 + res2)
        }

        return(1/E(Regressor, Afct, cent = cent, clip = clip, D1 = ErrorL2deriv))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "MultivariateDistribution", 
                                        neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, z.comp, clip, cent, stand, trafo){
        if(z.comp){
            wfct <- function(xx, u, clip, cent, stand){
                Y <- stand %*% (xx*u - cent)
                h.vct <- as.vector(sqrt(sum(Y^2)))
                ind2 <- (h.vct < clip/2)
                h.vct <- ind2*clip/2 + (1-ind2)*h.vct
                ind1 <- (h.vct < clip)

                return(ind1 + (1-ind1)*clip/h.vct)            
            }
            Afctu1 <- function(u, xx, clip, cent, stand){
                wfct(xx = xx, u = u, clip = clip, cent = cent, stand = stand)
            }
            Afctu2 <- function(u, xx, clip, cent, stand){
                u*wfct(xx = xx, u = u, clip = clip, cent = cent, stand = stand)
            }
            Afctu3 <- function(u, xx, clip, cent, stand){
                u^2*wfct(xx = xx, u = u, clip = clip, cent = cent, stand = stand)
            }
            Afct <- function(x, clip, cent, stand, D1){
                int1 <- E(D1, Afctu1, xx = x, clip = clip, cent = cent, stand = stand) 
                int2 <- E(D1, Afctu2, xx = x, clip = clip, cent = cent, stand = stand) 
                int3 <- E(D1, Afctu3, xx = x, clip = clip, cent = cent, stand = stand) 

                return((x %*% t(x))*int3 - (cent %*% t(x))*int2 
                        - (x %*% t(cent))*int2 + (cent %*% t(cent))*int1)
            }
            res <- E(Regressor, Afct, clip = clip, cent = cent, stand = stand,
                     D1 = ErrorL2deriv)
        }else{
            Afct <- function(x, clip, stand, D1){
                v <- t(x) %*% stand
                v <- as.vector(sqrt(v %*% t(v)))
                c0 <- clip/v

                return((x %*% t(x))*(m2df(D1, c0) - m2df(D1, -c0) - c0*(m1df(D1, -c0) + m1df(D1, c0))))
            }

            res <- E(Regressor, Afct, clip = clip, stand = stand, D1 = ErrorL2deriv)
        }

        return(trafo %*% solve(res))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution", 
                                        neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, z.comp, clip, cent, stand, trafo){
        Afct <- function(x, cent, clip, D1){
            if(x == 0) return(0)
            c1 <- cent(x) - clip/abs(x)
            c2 <- cent(x) + clip/abs(x)

            return(x^2*(m2df(D1, c2) - m2df(D1, c1) + c1*m1df(D1, c1) - c2*m1df(D1, c2)))
        }
        return(as.matrix(trafo/E(Regressor, Afct, cent = cent, clip = clip, D1 = ErrorL2deriv)))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "MultivariateDistribution", 
                                        neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, z.comp, clip, cent, stand, trafo){
        if(z.comp){
            wfct <- function(xx, u, clip, cent, stand){
                Y <- stand %*% (xx*u - cent(xx))
                h.vct <- as.vector(sqrt(sum(Y^2)))
                ind2 <- (h.vct < clip/2)
                h.vct <- ind2*clip/2 + (1-ind2)*h.vct
                ind1 <- (h.vct < clip)

                return(ind1 + (1-ind1)*clip/h.vct)            
            }
            Afctu <- function(u, xx, clip, cent, stand){
                (u - cent(xx))^2*wfct(xx = xx, u = u, clip = clip, cent = cent, stand = stand)
            }
            Afct <- function(x, clip, cent, stand, D1){
                return((x %*% t(x))*E(D1, Afctu, xx = x, clip = clip, 
                                      cent = cent, stand = stand))
            }
            res <- E(Regressor, Afct, clip = clip, cent = cent, stand = stand,
                     D1 = ErrorL2deriv)
        }else{
            Afct <- function(x, clip, stand, D1){
                v <- t(x) %*% stand
                v <- as.vector(sqrt(v %*% t(v)))
                c0 <- clip/v

                return((x %*% t(x))*(m2df(D1, c0) - m2df(D1, -c0) - c0*(m1df(D1, -c0) + m1df(D1, c0))))
            }

            res <- E(Regressor, Afct, clip = clip, stand = stand, D1 = ErrorL2deriv)
        }

        return(trafo %*% solve(res))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "Distribution", 
                                        neighbor = "Av2CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, z.comp, clip, cent, stand, trafo){
        c1 <- cent - clip
        c2 <- cent + clip
        return(1/(m2df(ErrorL2deriv, c2) - m2df(ErrorL2deriv, c1)
                + c1*m1df(ErrorL2deriv, c1) - c2*m1df(ErrorL2deriv, c2)))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "UnivariateDistribution", 
                                        neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, z.comp, clip, cent, stand, trafo){
        Afct <- function(x, cent, clip, D1){
            if(x == 0) return(0)
            gx <- cent(x)/abs(x)
            cx <- (cent(x) + clip)/abs(x)

            return(x^2*(m2df(D1, cx) - m2df(D1, gx) + gx*m1df(D1, gx) - cx*m1df(D1, cx)))
        }
        return(as.matrix(trafo/E(Regressor, Afct, cent = cent, clip = clip, D1 = ErrorL2deriv)))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                        Regressor = "MultivariateDistribution", 
                                        neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, z.comp, clip, cent, stand, trafo){
        Afct <- function(x, cent, clip, stand, D1){
            k <- length(x)
            if(all(x == numeric(k))) return(matrix(0, nrow = k, ncol = k))

            v <- as.vector(sqrt(sum((stand %*% x)^2)))
            gx <- cent(x)/v
            cx <- (cent(x) + clip)/v

            return((x %*% t(x))*(m2df(D1, cx) - m2df(D1, gx) + gx*m1df(D1, gx) - cx*m1df(D1, cx)))
        }
        return(trafo %*% solve(E(Regressor, Afct, cent = cent, clip = clip, 
                                 stand = stand, D1 = ErrorL2deriv)))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "RealRandVariable",
                                        Regressor = "Distribution", 
                                        neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, ErrorDistr, A.comp,
             stand, clip, cent, trafo){
        k <- length(ErrorL2deriv)
        integrandAu1 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            L1 <- as.vector(evalRandVar(ErrorL2deriv, u))
            L <- c(xx*L1[1], L1[2:k0])
            Y <- as.vector(stand %*% (L - cent))
            h.vct <- sqrt(sum(Y^2))
            ind2 <- (h.vct < clip/2)
            h.vct <- ind2*clip/2 + (1-ind2)*h.vct
            ind1 <- (h.vct < clip)

            return(ind1 + (1-ind1)*clip/h.vct)
        }
        integrandAu111 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            return(as.vector(ErrorL2deriv@Map[[1]](u))^2*integrandAu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandAu112 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            return(as.vector(ErrorL2deriv@Map[[1]](u))*integrandAu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandAu113 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            return(integrandAu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandAx11 <- function(x, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0, p){
            int1 <- E(object = ErrorDistr, fun = integrandAu111, xx = x, ErrorL2deriv = ErrorL2deriv, 
                      stand = stand, cent = cent, clip = clip, k0 = k0)
            int2 <- E(object = ErrorDistr, fun = integrandAu112, xx = x, ErrorL2deriv = ErrorL2deriv, 
                      stand = stand, cent = cent, clip = clip, k0 = k0)
            int3 <- E(object = ErrorDistr, fun = integrandAu113, xx = x, ErrorL2deriv = ErrorL2deriv, 
                      stand = stand, cent = cent, clip = clip, k0 = k0)
            return((x %*% t(x))*int1 - (cent[1:p] %*% t(x))*int2 
                   -(x %*% t(cent[1:p]))*int2 + (cent[1:p] %*% t(cent[1:p]))*int3)
        }
        integrandAu1j <- function(u, xx, Y.j, z.j, ErrorL2deriv, stand, cent, clip, k0, p){
            L2deriv1 <- (xx * as.vector(ErrorL2deriv@Map[[1]](u)) - cent[1:p]) * (Y.j(u) - z.j)
            return(L2deriv1*integrandAu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandAx1j <- function(x, Y.j, z.j, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0, p){
            E(object = ErrorDistr, fun = integrandAu1j, xx = x, Y.j = Y.j, z.j = z.j,
              ErrorL2deriv = ErrorL2deriv, stand = stand, cent = cent, clip = clip, k0 = k0, p = p)
        }
        integrandAuij <- function(u, xx, Y.i, Y.j, z.i, z.j, ErrorL2deriv, stand, cent, clip, k0){
            return((Y.i(u)-z.i) * (Y.j(u)-z.j)*integrandAu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                                    stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandAxij <- function(x, Y.i, Y.j, z.i, z.j, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandAuij, xx = x, Y.i = Y.i, Y.j = Y.j, z.i = z.i, 
              z.j = z.j, ErrorL2deriv = ErrorL2deriv, stand = stand, cent = cent, clip = clip, k0 = k0)
        }

        p <- dimension(img(Regressor))
        nrvalues <- p + k - 1
        res <- matrix(0, ncol=nrvalues, nrow=nrvalues)
        if(A.comp[1,1])
            res[1:p, 1:p] <- E(object = Regressor, fun = integrandAx11, ErrorL2deriv = ErrorL2deriv,
                               ErrorDistr = ErrorDistr, stand = stand, cent = cent, clip = clip, 
                               k0 = k, p = p)
        for(j in 2:k){
            if(A.comp[1,j])
                res[1:p, (p+j-1)] <- E(object = Regressor, fun = integrandAx1j, Y.j = ErrorL2deriv@Map[[j]], 
                                       z.j = cent[j], ErrorL2deriv = ErrorL2deriv, ErrorDistr = ErrorDistr, 
                                       stand = stand, cent = cent, clip = clip, k0 = k, p = p)
        } 
        for(i in 2:k)
            for(j in i:k){
                if(A.comp[i,j])
                    res[(p+i-1), (p+j-1)] <- E(object = Regressor, fun = integrandAxij, Y.i = ErrorL2deriv@Map[[i]], 
                                         Y.j = ErrorL2deriv@Map[[j]], z.i = cent[p+i-1], z.j = cent[p+j-1], 
                                         ErrorL2deriv = ErrorL2deriv, ErrorDistr = ErrorDistr, stand = stand, 
                                         cent = cent, clip = clip, k0 = k)
            }
        res[col(res) < row(res)] <- res[col(res) > row(res)]

        return(trafo %*% solve(res))
    })
setMethod("getInfStandRegTS", signature(ErrorL2deriv = "RealRandVariable",
                                        Regressor = "Distribution", 
                                        neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, ErrorDistr, A.comp,
             stand, clip, cent, trafo){
        k <- length(ErrorL2deriv)
        integrandAu1 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            L1 <- as.vector(evalRandVar(ErrorL2deriv, u))
            z1 <- cent(xx)
            L <- c(xx*L1[1]-z1[1], L1[2:k0]-z1[2:k0])
            Y <- as.vector(stand %*% L)
            h.vct <- sqrt(sum(Y^2))
            ind2 <- (h.vct < clip/2)
            h.vct <- ind2*clip/2 + (1-ind2)*h.vct
            ind1 <- (h.vct < clip)

            return(ind1 + (1-ind1)*clip/h.vct)
        }
        integrandAu11 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            L2deriv1 <- (as.vector(ErrorL2deriv@Map[[1]](u)) - cent(xx)[1])^2
            return(L2deriv1*integrandAu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandAx11 <- function(x, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            (x %*% t(x))*E(object = ErrorDistr, fun = integrandAu11, xx = x, ErrorL2deriv = ErrorL2deriv, 
              stand = stand, cent = cent, clip = clip, k0 = k0)
        }
        integrandAu1j <- function(u, xx, Y.j, j, ErrorL2deriv, stand, cent, clip, k0){
            z <- cent(xx)
            L2deriv1 <- (as.vector(ErrorL2deriv@Map[[1]](u)) - z[1])*(Y.j(u) - z[j])
            return(L2deriv1*integrandAu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandAx1j <- function(x, Y.j, j, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            x * E(object = ErrorDistr, fun = integrandAu1j, xx = x, Y.j = Y.j, j = j,
                  ErrorL2deriv = ErrorL2deriv, stand = stand, cent = cent, clip = clip, k0 = k0)
        }
        integrandAuij <- function(u, xx, Y.i, Y.j, indi, indj, ErrorL2deriv, stand, cent, clip, k0){
            z <- cent(xx)
            return((Y.i(u)-z[indi]) * (Y.j(u)-z[indj])*integrandAu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                                    stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandAxij <- function(x, Y.i, Y.j, indi, indj, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandAuij, xx = x, Y.i = Y.i, Y.j = Y.j, indi = indi, 
              indj = indj, ErrorL2deriv = ErrorL2deriv, stand = stand, cent = cent, clip = clip, k0 = k0)
        }

        p <- dimension(img(Regressor))
        nrvalues <- p + k - 1
        res <- matrix(0, ncol=nrvalues, nrow=nrvalues)
        if(A.comp[1,1])
            res[1:p, 1:p] <- E(object = Regressor, fun = integrandAx11, ErrorL2deriv = ErrorL2deriv,
                               ErrorDistr = ErrorDistr, stand = stand, cent = cent, clip = clip, 
                               k0 = k)
        for(j in 2:k){
            if(A.comp[1,j])
                res[1:p, (p+j-1)] <- E(object = Regressor, fun = integrandAx1j, Y.j = ErrorL2deriv@Map[[j]], 
                                       j = j, ErrorL2deriv = ErrorL2deriv, ErrorDistr = ErrorDistr, 
                                       stand = stand, cent = cent, clip = clip, k0 = k)
        } 
        for(i in 2:k)
            for(j in i:k){
                if(A.comp[i,j])
                    res[(p+i-1), (p+j-1)] <- E(object = Regressor, fun = integrandAxij, Y.i = ErrorL2deriv@Map[[i]], 
                                         Y.j = ErrorL2deriv@Map[[j]], indi = i, indj = j, ErrorL2deriv = ErrorL2deriv, 
                                         ErrorDistr = ErrorDistr, stand = stand, cent = cent, clip = clip, k0 = k)
            }
        res[col(res) < row(res)] <- res[col(res) > row(res)]

        return(trafo %*% solve(res))
    })
