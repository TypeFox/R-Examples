###############################################################################
## centering constant for asymptotic MSE and asymptotic Hampel risk 
## in case of (unconditional) contamination neighborhoods (* = c)
###############################################################################
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "UnivariateDistribution",
                                       neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, stand, z.comp){
        if(!z.comp) return(numeric(1))
        
        z.fct <- function(x, z, c0, D1){
            if(x == 0) return(0)
            c0.x <- c0/abs(x)
            z.x <- z/x
            
            return(x*((c0.x*(1 - p(D1)(z.x-c0.x) - p(D1)(z.x+c0.x)) 
                    + m1df(D1, z.x+c0.x) - m1df(D1, z.x-c0.x))/
                   (p(D1)(z.x+c0.x) - p(D1)(z.x-c0.x))))
        }
        return(E(Regressor, z.fct, z = cent, c0 = clip, D1 = ErrorL2deriv))        
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "UnivariateDistribution",
                                       neighbor = "TotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, z.comp){
        if(!z.comp) return(-clip/2)
        
        g.fct <- function(z, c0, D1, K){
            gu.fct <- function(x, z, c0, D1){
                if(x == 0) return(0)
                c1 <- (z + c0)/x
                c2 <- z/x
                res1 <- (- max(x,0)*(m1df(D1, c1) + c1*(1-p(D1)(c1)))
                         + min(x, 0)*(m1df(D1, c1) - c1*p(D1)(c1)))
                res2 <- (max(x,0)*(c2*p(D1)(c2) - m1df(D1, c2)) 
                         + min(x,0)*(c2*(1-p(D1)(c2)) + m1df(D1, c2)))

                return(res2 - res1)
            }
            return(E(K, gu.fct, z = z, c0 = c0, D1 = D1))
        }
        lower <- q(ErrorL2deriv)(distr::TruncQuantile)
        upper <- q(ErrorL2deriv)(1-distr::TruncQuantile)

        return(uniroot(g.fct, lower = lower, upper = upper, tol = tol.z, 
                    c0 = clip, D1 = D1, K = Regressor)$root)        
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "MultivariateDistribution",
                                       neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, stand, z.comp){
        if(!z.comp) return(numeric(ncol(stand)))

        integrandzu1 <- function(u, xx, stand, cent, clip){
            Y <- stand %*% (xx*u - cent)
            h.vct <- as.vector(sqrt(sum(Y^2)))
            ind2 <- (h.vct < clip/2)
            h.vct <- ind2*clip/2 + (1-ind2)*h.vct
            ind1 <- (h.vct < clip)

            return(ind1 + (1-ind1)*clip/h.vct)
        }
        integrandzx1 <- function(x, ErrorL2deriv, stand, cent, clip){
            E(object = ErrorL2deriv, fun = integrandzu1, xx = x,  
              stand = stand, cent = cent, clip = clip)
        }
        integrandzu2 <- function(u, xx, stand, cent, clip){
            return(u*integrandzu1(u = u, xx = xx, stand = stand, cent = cent, clip = clip))
        }
        integrandzx2 <- function(x, ErrorL2deriv, stand, cent, clip){
            as.vector(x)*E(object = ErrorL2deriv, fun = integrandzu2, xx = x,
                           stand = stand, cent = cent, clip = clip)
        }
        
        res1 <- E(object = Regressor, fun = integrandzx1, ErrorL2deriv = ErrorL2deriv, 
                  stand = stand, cent = cent, clip = clip)
        res2 <- E(object = Regressor, fun = integrandzx2, ErrorL2deriv = ErrorL2deriv, 
                  stand = stand, cent = cent, clip = clip)

        return(res2/res1)
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "numeric",
                                       neighbor = "CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, z.comp){
        if(!z.comp) return(-clip)
        
        if(x > 0){
            g.fct <- function(z, c0, D1, x){
                z*p(D1)(z/x) - x*(m1df(D1, z/x) - m1df(D1, b/x)) + b*(1-p(D1)(b/x))
            }
        }else{
            g.fct <- function(z, c0, D1, x){
                z*(1-p(D1)(z/x)) + x*(m1df(D1, z/x) - m1df(D1, b/x)) + b*p(D1)(b/x)
            }
        }
        lower <- q(ErrorL2deriv)(distr::TruncQuantile)
        upper <- q(ErrorL2deriv)(1-distr::TruncQuantile)

        return(uniroot(g.fct, lower = lower, upper = upper, tol = tol.z, 
                    c0 = clip, D1 = D1, x = Regressor)$root)        
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "UnivariateDistribution",
                                       neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, stand, z.comp, x.vec){
        if(!z.comp){ 
            z.fct <- function(x){numeric(1)}
            z.vec <- numeric(length(x.vec))
            return(list(z.fct = z.fct, z.vec = z.vec))
        }
        
        zfun <- function(x, z0, c0, D1){
            if(x == 0) return(0)
            c0.x <- c0/abs(x)
            z.x <- z0(x)
            
            return(x*((c0.x*(1 - p(D1)(z.x-c0.x) - p(D1)(z.x+c0.x)) 
                    + m1df(D1, z.x+c0.x) - m1df(D1, z.x-c0.x))/
                   (p(D1)(z.x+c0.x) - p(D1)(z.x-c0.x))))
        }
        z.vec <- sapply(x.vec, zfun, z0 = cent, c0 = clip, D1 = ErrorL2deriv)

        if(is(Regressor, "DiscreteDistribution")){
            z.fct <- function(x){ 
                ind <- (round(x, 8) == round(x.vec, 8))
                if(any(ind))
                    return(z.vec[ind])
                else
                    return(NA)
            }
        }else{
            if(is.finite(q(Regressor)(0)))
                yleft <- NA
            else
                yleft <- z.vec[1]
            if(is.finite(q(Regressor)(1)))
                yright <- NA
            else
                yright <- z.vec[length(z.vec)]
            z.fct <- approxfun(x.vec, z.vec, yleft = yleft, yright = yright)
        }

        return(list(z.fct = z.fct, z.vec = z.vec))
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "MultivariateDistribution",
                                       neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, stand, z.comp, x.vec){
        if(!z.comp){ 
            z.fct <- function(x){numeric(1)}
            z.vec <- numeric(nrow(as.matrix(x.vec)))
            return(list(z.fct = z.fct, z.vec = z.vec))
        }

        integrandzu1 <- function(u, xx, stand, cent, clip){
            Y <- stand %*% (xx*u - cent(xx))
            h.vct <- as.vector(sqrt(sum(Y^2)))
            ind2 <- (h.vct < clip/2)
            h.vct <- ind2*clip/2 + (1-ind2)*h.vct
            ind1 <- (h.vct < clip)

            return(ind1 + (1-ind1)*clip/h.vct)
        }
        integrandzx1 <- function(x, ErrorL2deriv, stand, cent, clip){
            E(object = ErrorL2deriv, fun = integrandzu1, xx = x,  
              stand = stand, cent = cent, clip = clip)
        }
        integrandzu2 <- function(u, xx, stand, cent, clip){
            return(u*integrandzu1(u = u, xx = xx, stand = stand, cent = cent, clip = clip))
        }
        integrandzx2 <- function(x, ErrorL2deriv, stand, cent, clip){
            E(object = ErrorL2deriv, fun = integrandzu2, xx = x,
              stand = stand, cent = cent, clip = clip)
        }
        
        res1 <- as.vector(apply(x.vec, 1, integrandzx1, ErrorL2deriv = ErrorL2deriv, 
                      stand = stand, cent = cent, clip = clip))
        res2 <- as.vector(apply(x.vec, 1, integrandzx2, ErrorL2deriv = ErrorL2deriv,
                      stand = stand, cent = cent, clip = clip))
        z.vec <- res2/res1
        k <- dimension(img(Regressor))
        if(is(Regressor, "DiscreteMVDistribution")){
            z.fct <- function(x){ 
                ind <- colSums(apply(round(x.vec, 8), 1, "==", round(x, 8))) == k
                if(any(ind))
                    return(z.vec[ind])
                else
                    return(NA)
            }
        }else{
            stop("not yet implemented")
        }

        return(list(z.fct = z.fct, z.vec = z.vec))
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "Distribution",
                                       neighbor = "Av2CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, stand, z.comp, tol.z){
        if(!z.comp) return(0)

        z.fct <- function(z, c0, D1){
            return(c0 + (z-c0)*p(D1)(z-c0) - (z+c0)*p(D1)(z+c0) + m1df(D1, z+c0) - m1df(D1, z-c0))
        }
        lower <- q(ErrorL2deriv)(distr::TruncQuantile)
        upper <- q(ErrorL2deriv)(1-distr::TruncQuantile)

        return(uniroot(z.fct, lower = lower, upper = upper, tol = tol.z, 
                    c0=clip, D1=ErrorL2deriv)$root)
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "UnivariateDistribution",
                                       neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, stand, z.comp, x.vec, tol.z){
        if(!z.comp){ 
            z.fct <- function(x){-b/2}
            body(z.fct) <- substitute({-b/2}, list(b = clip))
            z.vec <- numeric(length(x.vec)) - clip/2
            return(list(z.fct = z.fct, z.vec = z.vec))
        }
        
        g.fct <- function(g, c0, xx, D1){
            gx <- g/abs(xx)
            cx <- c0/abs(xx)
            return(gx*p(D1)(gx) + (gx+cx)*(1-p(D1)(gx+cx)) - m1df(D1, gx) + m1df(D1, gx+cx))
        }
        zfun <- function(x, z0, c0, D1, tol.z){
            if(x == 0) return(0)
            
            lower <- q(D1)(distr::TruncQuantile)
            upper <- q(D1)(1-distr::TruncQuantile)

            return(uniroot(g.fct, lower = lower, upper = upper, tol = tol.z, 
                        c0 = c0, xx = x, D1 = D1)$root)
        }
        z.vec <- sapply(x.vec, zfun, z0 = cent, c0 = clip, D1 = ErrorL2deriv, tol.z = tol.z)

        if(is(Regressor, "DiscreteDistribution")){
            z.fct <- function(x){ 
                ind <- (round(x, 8) == round(x.vec, 8))
                if(any(ind))
                    return(z.vec[ind])
                else
                    return(NA)
            }
        }else{
            if(is.finite(q(Regressor)(0)))
                yleft <- NA
            else
                yleft <- z.vec[1]
            if(is.finite(q(Regressor)(1)))
                yright <- NA
            else
                yright <- z.vec[length(z.vec)]
            z.fct <- approxfun(x.vec, z.vec, yleft = yleft, yright = yright)
        }

        return(list(z.fct = z.fct, z.vec = z.vec))
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "MultivariateDistribution",
                                       neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, clip, cent, stand, z.comp, x.vec, tol.z){
        if(!z.comp){ 
            z.fct <- function(x){ -b/2 }
            body(z.fct) <- substitute({ -b/2 }, list(b = clip))
            z.vec <- numeric(nrow(as.matrix(x.vec))) - clip/2
            return(list(z.fct = z.fct, z.vec = z.vec))
        }

        g.fct <- function(g, c0, A0, xx, D1){
            v <- as.vector(sqrt(sum((A0 %*% xx)^2)))
            gx <- g/v
            cx <- c0/v
            return(gx*p(D1)(gx) + (gx+cx)*(1-p(D1)(gx+cx)) - m1df(D1, gx) + m1df(D1, gx+cx))
        }
        zfun <- function(x, z0, c0, A0, D1, tol.z){
            if(all(x == numeric(length(x)))) return(0)
            
            lower <- q(D1)(distr::TruncQuantile)
            upper <- q(D1)(1-distr::TruncQuantile)

            return(uniroot(g.fct, lower = lower, upper = upper, tol = tol.z, 
                        c0 = c0, A0 = A0, xx = x, D1 = D1)$root)
        }
        z.vec <- apply(x.vec, 1, zfun, z0 = cent, c0 = clip, A0 = stand, 
                        D1 = ErrorL2deriv, tol.z = tol.z)

        k <- dimension(img(Regressor))
        if(is(Regressor, "DiscreteMVDistribution")){
            z.fct <- function(x){ 
                ind <- colSums(apply(round(x.vec, 8), 1, "==", round(x, 8))) == k
                if(any(ind))
                    return(z.vec[ind])
                else
                    return(NA)
            }
        }else{
            stop("not yet implemented")
        }

        return(list(z.fct = z.fct, z.vec = z.vec))
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "RealRandVariable",
                                       Regressor = "Distribution",
                                       neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, ErrorDistr, stand, cent, clip, z.comp){
        integrandzu1 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            L1 <- as.vector(evalRandVar(ErrorL2deriv, u))
            L <- c(as.vector(xx)*L1[1], L1[2:k0])
            Y <- as.vector(stand %*% (L - cent))
            h.vct <- sqrt(sum(Y^2))
            ind2 <- (h.vct < clip/2)
            h.vct <- ind2*clip/2 + (1-ind2)*h.vct
            ind1 <- (h.vct < clip)

            return(ind1 + (1-ind1)*clip/h.vct)
        }
        integrandzx1 <- function(x, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandzu1, xx = x, ErrorL2deriv = ErrorL2deriv, 
              stand = stand, cent = cent, clip = clip, k0 = k0)
        }
        integrandzu2 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            L2deriv1 <- as.vector(ErrorL2deriv@Map[[1]](u))
            return(L2deriv1*integrandzu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandzx2 <- function(x, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            as.vector(x)*E(object = ErrorDistr, fun = integrandzu2, xx = x, ErrorL2deriv = ErrorL2deriv, 
              stand = stand, cent = cent, clip = clip, k0 = k0)
        }
        integrandzu3 <- function(u, xx, Y.i, ErrorL2deriv, stand, cent, clip, k0){
            return(Y.i(u)*integrandzu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandzx3 <- function(x, Y.i, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandzu3, xx = x, Y.i = Y.i, 
              ErrorL2deriv = ErrorL2deriv, stand = stand, cent = cent, clip = clip, k0 = k0)
        }
        
        k <- length(ErrorL2deriv)
        res1 <- E(object = Regressor, fun = integrandzx1, ErrorL2deriv = ErrorL2deriv, 
                  ErrorDistr = ErrorDistr, stand = stand, cent = cent, clip = clip, k0 = k)
        p <- dimension(img(Regressor))
        nrvalues <- p + k - 1
        res2 <- numeric(nrvalues)
        if(z.comp[1]){
            res2[1:p] <- E(object = Regressor, fun = integrandzx2, ErrorL2deriv = ErrorL2deriv, 
                           ErrorDistr = ErrorDistr, stand = stand, cent = cent, clip = clip, k0 = k)
        }
        for(i in 2:k){
            if(z.comp[i]){
                res2[p+i-1] <- E(object = Regressor, fun = integrandzx3, Y.i = ErrorL2deriv@Map[[i]], 
                                 ErrorL2deriv = ErrorL2deriv, ErrorDistr = ErrorDistr, stand = stand, 
                                 cent = cent, clip = clip, k0 = k)
            }
        }

        return(res2/res1)
    })
setMethod("getInfCentRegTS", signature(ErrorL2deriv = "RealRandVariable",
                                       Regressor = "Distribution",
                                       neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, neighbor, ErrorDistr, stand, cent, clip, z.comp, x.vec){
        integrandzu1 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
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
        integrandzx1 <- function(x, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandzu1, xx = x, ErrorL2deriv = ErrorL2deriv, 
              stand = stand, cent = cent, clip = clip, k0 = k0)
        }
        integrandzu2 <- function(u, xx, ErrorL2deriv, stand, cent, clip, k0){
            L2deriv1 <- as.vector(ErrorL2deriv@Map[[1]](u))
            return(L2deriv1*integrandzu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandzx2 <- function(x, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandzu2, xx = x, ErrorL2deriv = ErrorL2deriv, 
              stand = stand, cent = cent, clip = clip, k0 = k0)
        }
        integrandzu3 <- function(u, xx, Y.i, ErrorL2deriv, stand, cent, clip, k0){
            return(Y.i(u)*integrandzu1(u = u, xx = xx, ErrorL2deriv = ErrorL2deriv, 
                                       stand = stand, cent = cent, clip = clip, k0 = k0))
        }
        integrandzx3 <- function(x, Y.i, ErrorL2deriv, ErrorDistr, stand, cent, clip, k0){
            E(object = ErrorDistr, fun = integrandzu3, xx = x, Y.i = Y.i, 
              ErrorL2deriv = ErrorL2deriv, stand = stand, cent = cent, clip = clip, k0 = k0)
        }
        
        k <- length(ErrorL2deriv)
        if(is(Regressor, "UnivariateDistribution")){
            res1 <- sapply(x.vec, integrandzx1, ErrorL2deriv = ErrorL2deriv, 
                           ErrorDistr = ErrorDistr, stand = stand, cent = cent, 
                           clip = clip, k0 = k)
            res2 <- matrix(0, ncol = length(x.vec), nrow = k)
            if(z.comp[1])
                res2[1,] <- sapply(x.vec, integrandzx2, ErrorL2deriv = ErrorL2deriv, 
                                   ErrorDistr = ErrorDistr, stand = stand, cent = cent, 
                                   clip = clip, k0 = k)
            for(i in 2:k)
                if(z.comp[i])
                    res2[i,] <- sapply(x.vec, integrandzx3, Y.i = ErrorL2deriv@Map[[i]], 
                                       ErrorL2deriv = ErrorL2deriv, ErrorDistr = ErrorDistr, 
                                       stand = stand, cent = cent, clip = clip, k0 = k)
            z.vec <- t(res2)/res1

            if(is(Regressor, "DiscreteDistribution")){
                z.fct <- function(x){ 
                    ind <- (round(x, 8) == round(x.vec, 8))
                    if(any(ind))
                        return(z.vec[ind,])
                    else
                        return(NA)
                }
            }else{
                z.fun <- vector("list", length = k)
                for(i in 1:k)
                    z.fun[[i]] <- approxfun(x.vec, z.vec[,i], rule = 2)

                z.fct <- function(x){
                    z <- numeric(k)
                    for(i in 1:k) z[i] <- z.fun[[i]](x)
                    return(z)
                }
            }
        }else{
            res1 <- apply(x.vec, 1, integrandzx1, ErrorL2deriv = ErrorL2deriv, 
                           ErrorDistr = ErrorDistr, stand = stand, cent = cent, 
                           clip = clip, k0 = k)
            res2 <- matrix(0, ncol = nrow(x.vec), nrow = k)
            if(z.comp[1]){
                res2[1,] <- apply(x.vec, 1, integrandzx2, ErrorL2deriv = ErrorL2deriv, 
                                  ErrorDistr = ErrorDistr, stand = stand, cent = cent, 
                                  clip = clip, k0 = k)
            }
            for(i in 2:k){
                if(z.comp[i]){
                    res2[i,] <- apply(x.vec, 1, integrandzx3, Y.i = ErrorL2deriv@Map[[i]], 
                                      ErrorL2deriv = ErrorL2deriv, ErrorDistr = ErrorDistr, 
                                      stand = stand, cent = cent, clip = clip, k0 = k)
                }
            }
            z.vec <- t(res2)/res1
            
            z.fct <- function(x){ 
                ind <- colSums(apply(round(x.vec, 8), 1, "==", round(x, 8))) == k
                if(any(ind))
                    return(z.vec[ind,])
                else
                    return(NA)
            }
        }
        return(list(z.fct = z.fct, z.vec = z.vec))
    })
