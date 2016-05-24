###############################################################################
## asymptotic MSE
###############################################################################
setMethod("getAsRiskRegTS", signature(risk = "asMSE",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "Distribution",
                                      neighbor = "Neighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, clip, cent, stand, trafo){
        if(!is.finite(neighbor@radius))
            return(list(asMSE = Inf))
        else
            return(list(asMSE = sum(diag(stand %*% t(trafo)))))
    })
setMethod("getAsRiskRegTS", signature(risk = "asMSE",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "Distribution",
                                      neighbor = "Av2CondContNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, clip, cent, stand, trafo){
        if(!is.finite(neighbor@radius))
            return(list(asMSE = Inf))
        else{
            K.inv <- solve(E(Regressor, fun = function(x){ x %*% t(x) }))
            return(list(asMSE = stand * sum(diag(t(trafo) %*% K.inv))))
        }
    })
setMethod("getAsRiskRegTS", signature(risk = "asMSE",
                                      ErrorL2deriv = "EuclRandVariable",
                                      Regressor = "Distribution",
                                      neighbor = "Neighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, clip, cent, stand, trafo){
        if(!is.finite(neighbor@radius))
            return(list(asMSE = Inf))
        else
            return(list(asMSE = sum(diag(stand %*% t(trafo)))))
    })

###############################################################################
## standardized asymptotic bias
###############################################################################
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "ContNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, 
             ErrorL2derivDistrSymm, trafo, maxiter, tol){
        if(is(ErrorL2derivDistrSymm, "SphericalSymmetry")) 
            z.comp <- !(ErrorL2derivDistrSymm@SymmCenter == 0)
        else
            z.comp <- TRUE
        
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        if(z.comp){
            abs.fct <- function(x, ErrorL2deriv, cent){ 
                abs(x)*E(ErrorL2deriv, function(u, xx, cent){abs(u - cent/x)}, xx = x, cent = cent)
            }
            bmin.fct <- function(cent, ErrorL2deriv, Regressor){
                E(Regressor, abs.fct, ErrorL2deriv = ErrorL2deriv, cent = cent)
            }
            erg <- optimize(bmin.fct, lower = -upper, upper = upper, 
                            tol = tol, ErrorL2deriv = ErrorL2deriv, Regressor = Regressor)
            b <- as.vector(trafo)/erg$objective
        }else{
            b <- zi*as.vector(trafo)/(-2*m1df(ErrorL2deriv, 0)*E(Regressor, function(x){abs(x)}))
        }

        return(list(asBias = b))        
    })
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "Av1CondContNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, 
             ErrorL2derivDistrSymm, trafo, maxiter, tol){
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        z <- q(ErrorL2deriv)(0.5)
        Eu <- E(ErrorL2deriv, function(x, z){abs(x - z)}, z = z)
        Ex <- E(Regressor, function(x){abs(x)})
        b <- zi*as.vector(trafo)/(Ex*Eu)

        return(list(asBias = b))        
    })
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "Av1CondTotalVarNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, 
             ErrorL2derivDistrSymm, trafo, maxiter, tol){
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        Ex <- E(Regressor, fun = function(x){abs(x)})

        return(zi*as.vector(trafo)/(-m1df(ErrorL2deriv, 0)*Ex))
    })
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "MultivariateDistribution",
                                      neighbor = "ContNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, 
             ErrorL2derivDistrSymm, trafo, maxiter, tol){
        if(is(ErrorL2derivDistrSymm, "SphericalSymmetry")) 
            z.comp <- !(ErrorL2derivDistrSymm@SymmCenter == 0)
        else
            z.comp <- TRUE

        if(z.comp){
            abs.fctu <- function(u, xx, A, a0){ 
                v <- t(A %*% xx * u - a0)
                return(as.vector(sqrt(v %*% t(v))))
            }
            abs.fctx <- function(x, ErrorL2deriv, A, a0){
                E(ErrorL2deriv, abs.fctu, xx = x, A = A, a0 = a0)
            }

            bmin.fct <- function(param, ErrorL2deriv, Regressor, trafo){
                p <- nrow(trafo)
                k <- ncol(trafo)
                A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
                a <- param[(p*k+1):length(param)]
    
                return(E(Regressor, abs.fctx, ErrorL2deriv = ErrorL2deriv, A = A, a0 = a)/sum(diag(A %*% t(trafo))))
            }
        
            erg <- optim(c(as.vector(trafo), numeric(nrow(trafo))), bmin.fct, method = "Nelder-Mead", 
                        control=list(reltol=tol, maxit=100*maxiter), Regressor = Regressor, 
                        ErrorL2deriv = ErrorL2deriv, trafo = trafo)
            b <- 1/erg$value
        }else{
            abs.fct <- function(x, A, a) {
                v <- t(A %*% x)
                return(as.vector(sqrt(v %*% t(v))))
            }
            bmin.fct <- function(param, Regressor, trafo) {
                A <- matrix(param, ncol = ncol(trafo), nrow = nrow(trafo))
                return(E(Regressor, abs.fct, A = A)/sum(diag(A %*% t(trafo))))
            }
            erg <- optim(as.vector(trafo), bmin.fct, method = "Nelder-Mead", 
                         control = list(reltol = tol, maxit = 100 * maxiter), 
                         Regressor = Regressor, trafo = trafo)
            b <- 1/(erg$value * E(ErrorL2deriv, function(x) { abs(x) }))
        }
    
        return(list(asBias = b))        
    })
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "MultivariateDistribution",
                                      neighbor = "Av1CondContNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, 
             ErrorL2derivDistrSymm, trafo, maxiter, tol){
        abs.fct <- function(x, A){ 
            v <- t(x) %*% t(A)
            return(as.vector(sqrt(v %*% t(v))))
        }

        bmin.fct <- function(param, Regressor, trafo){
            A <- matrix(param, ncol = ncol(trafo), nrow = nrow(trafo))

            return(E(Regressor, abs.fct, A = A)/sum(diag(A %*% t(trafo))))
        }
        
        erg <- optim(as.vector(trafo), bmin.fct, method = "Nelder-Mead", 
                     control = list(reltol = tol, maxit = 100*maxiter), 
                     Regressor = Regressor, trafo = trafo)
        z <- q(ErrorL2deriv)(0.5)
        b <- 1/(erg$value*E(ErrorL2deriv, function(x, z){abs(x - z)}, z = z))

        return(list(asBias = b))        
    })
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "Distribution",
                                      neighbor = "Av2CondContNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, 
             ErrorL2derivDistrSymm, trafo, maxiter, tol){
        K <- E(Regressor, fun = function(x){ x %*% t(x) })
        z <- q(ErrorL2deriv)(0.5)
        Eu <- E(ErrorL2deriv, function(x, z){abs(x - z)}, z = z)
        b <- sqrt(sum(diag(trafo %*% solve(K) %*% t(trafo))))/Eu
        
        return(list(asBias = b))        
    })
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "MultivariateDistribution",
                                      neighbor = "Av1CondTotalVarNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, 
             ErrorL2derivDistrSymm, trafo, maxiter, tol){
        abs.fct <- function(x, A){ 
            return(as.vector(sqrt(sum((A %*% x)^2))))
        }

        bmin.fct <- function(param, Regressor, trafo, Eu){
            A <- matrix(param, ncol = ncol(trafo), nrow = nrow(trafo))
            return(E(Regressor, abs.fct, A = A)/sum(diag(A %*% t(trafo))))
        }
        
        erg <- optim(as.vector(trafo), bmin.fct, method = "Nelder-Mead", 
                    control=list(reltol=tol, maxit=100*maxiter), Regressor = Regressor, trafo = trafo)

        return(1/(erg$value*(-m1df(ErrorL2deriv, 0))))
    })
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "RealRandVariable",
                                      Regressor = "Distribution",
                                      neighbor = "ContNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, ErrorDistr, trafo, 
             z.start, A.start, maxiter, tol){
        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        abs.fctu <- function(u, xx, ErrorL2deriv, A, a0, k){ 
            L1 <- as.vector(evalRandVar(ErrorL2deriv, u))
            L2deriv <- c(xx*L1[1], L1[2:k])
            Y <- as.vector(A %*% L2deriv - a0)

            return(sqrt(sum(Y^2)))
        }
        abs.fctx <- function(x, ErrorL2deriv, ErrorDistr, A, a0, k){ 
            E(object = ErrorDistr, fun = abs.fctu, xx = x, ErrorL2deriv = ErrorL2deriv, 
              A = A, a0 = a0, k = k)
        }

        bmin.fct <- function(param, ErrorL2deriv, Regressor, ErrorDistr, trafo){
            p <- nrow(trafo)
            k <- ncol(trafo)
            A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
            a <- param[(p*k+1):length(param)]

            return(E(object = Regressor, fun = abs.fctx, ErrorL2deriv = ErrorL2deriv,
                     ErrorDistr = ErrorDistr, A = A, a0 = a, k = k)/sum(diag(A %*% t(trafo))))
        }
        
        A.vec <- as.vector(A.start)
        erg <- optim(c(A.vec, as.vector(A.start %*% z.start)), bmin.fct, method = "Nelder-Mead", 
                    control=list(reltol=tol, maxit=100*maxiter), 
                    ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, ErrorDistr = ErrorDistr, 
                    trafo = trafo)
        b <- 1/erg$value

        return(list(asBias = b))        
    })
setMethod("getAsRiskRegTS", signature(risk = "asBias",
                                      ErrorL2deriv = "RealRandVariable",
                                      Regressor = "Distribution",
                                      neighbor = "Av1CondContNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, ErrorDistr, trafo, 
             z.start, A.start, maxiter, tol){
        if(is.null(z.start)) z.start <- function(x){numeric(1)}
        if(is.null(A.start)) A.start <- trafo
        
        stop("not yet implemented")

        return(list(asBias = b))        
    })
setMethod("getAsRiskRegTS", signature(risk = "asUnOvShoot",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "UncondNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, clip, cent, stand){
        if(identical(all.equal(neighbor@radius, 0), TRUE))
            return(list(asUnOvShoot = pnorm(-risk@width/sqrt(as.vector(stand)))))
        
        var.fct <- function(x, cent, clip, D1){
            g0 <- cent/x
            c0 <- (cent + clip)/x
            res1 <- max(0,x)^2*(g0^2*p(D1)(g0) + m2df(D1, c0) - m2df(D1, g0) 
                                + c0^2*(1-p(D1)(c0)))
            res2 <- min(0,x)^2*(c0^2*p(D1)(c0) + m2df(D1, g0) - m2df(D1, c0)
                                + g0^2*(1-p(D1)(g0)))
            return(res1 + res2)
        }
        s <- sqrt(E(Regressor, var.fct, cent = cent, clip = clip, D1 = ErrorL2deriv))

        return(list(asUnOvShoot = pnorm(-risk@width*s)))
    })
setMethod("getAsRiskRegTS", signature(risk = "asUnOvShoot",
                                      ErrorL2deriv = "UnivariateDistribution",
                                      Regressor = "UnivariateDistribution",
                                      neighbor = "CondNeighborhood"),
    function(risk, ErrorL2deriv, Regressor, neighbor, clip, cent, stand){
        if(is(Regressor, "AbscontDistribution")){
            xlower <- ifelse(is.finite(q(Regressor)(0)), q(Regressor)(0), q(Regressor)(distr::TruncQuantile))
            xupper <- ifelse(is.finite(q(Regressor)(1)), q(Regressor)(1), q(Regressor)(1 - distr::TruncQuantile))
            x.vec <- seq(from = xlower, to = xupper, by = 0.01)
        }else{
            if(is(Regressor, "DiscreteDistribution"))
                x.vec <- support(Regressor) 
            else
                x.vec <- unique(r(Regressor)(1e3))
        }
        if(identical(all.equal(max(neighbor@radiusCurve(x.vec)), 0), TRUE))
            return(list(asUnOvShoot = pnorm(-risk@width/sqrt(as.vector(stand)))))
        
        var.fct <- function(x, cent, clip, D1){
            g0 <- cent(x)/x
            c0 <- clip(x)/x
            res1 <- max(0,x)^2*(g0^2*p(D1)(g0) + m2df(D1, c0) - m2df(D1, g0) 
                                + c0^2*(1-p(D1)(c0)))
            res2 <- min(0,x)^2*(c0^2*p(D1)(c0) + m2df(D1, g0) - m2df(D1, c0)
                                + g0^2*(1-p(D1)(g0)))
            return(res1 + res2)
        }
        s <- sqrt(E(Regressor, var.fct, cent = cent, clip = clip, D1 = ErrorL2deriv))

        return(list(asUnOvShoot = pnorm(-risk@width*s)))
    })
