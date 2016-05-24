###############################################################################
## get minimum bias solutions
###############################################################################
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "UnivariateDistribution",
                                          risk = "asBias", 
                                          neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo, upper, maxiter, tol, warn){
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
            a <- as.vector(A %*% erg$minimum)
            
            if(is(ErrorL2deriv, "AbscontDistribution") && is(Regressor, "AbscontDistribution")){
                d <- 0
            }

            if(is(ErrorL2deriv, "AbscontDistribution") && is(Regressor, "DiscreteDistribution"))
                if(identical(all.equal(erg$minimum, 0, tolerance = sqrt(tol)), TRUE) && d(Regressor)(0) == 0)
                    d <- 0
                else{
                    d <- 0
                    cat("computation of 'd' not yet implemented\n")
                }
            
            info <- c("minimum asymptotic bias (lower case) solution")
            asCov <- b^2
            Risk <- list(asBias = b, asCov = asCov)
        }else{
            b <- zi*as.vector(trafo)/(-2*m1df(ErrorL2deriv, 0)*E(Regressor, function(x){abs(x)}))

            if(is(Regressor, "AbscontDistribution") || (is(Regressor, "DiscreteDistribution") && d(Regressor)(0) == 0)){
                if(is(ErrorL2deriv, "AbscontDistribution")){
                    d <- 0
                }else{
                    if(d(ErrorL2deriv)(0) == 0)
                        d <- 0
                    else{
                        d <- 0
                        cat("computation of 'd' not yet implemented\n")
                    }
                }
            }else{
                d <- 0
                cat("computation of 'd' not yet implemented\n")
            }
            
            a <- 0
            info <- c("minimum asymptotic bias (lower case) solution")
            asCov <- b^2
            Risk <- list(asBias = b, asCov = asCov)
        }

        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info))    
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "UnivariateDistribution",
                                          risk = "asBias", 
                                          neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo, upper, maxiter, tol, warn){
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        z <- q(ErrorL2deriv)(0.5)
        Eu <- E(ErrorL2deriv, function(x, z){abs(x - z)}, z = z)
        Ex <- E(Regressor, abs)
        b <- zi*as.vector(trafo)/(Ex*Eu)

        if(is(Regressor, "AbscontDistribution") | (is(Regressor, "DiscreteDistribution") & d(Regressor)(0) == 0)){
            if(is(ErrorL2deriv, "AbscontDistribution")){
                ws0 <- 0
                d <- 0
            }else{
                ws0 <- d(ErrorL2deriv)(z)
                d <- (2*p(ErrorL2deriv)(z) - ws0 - 1)/ws0                
            }
        }else{
            d <- 0
            ws0 <- 0
            cat("computation of 'd' not yet implemented\n")
        }

        info <- c("minimum asymptotic bias (lower case) solution")
        asCov <- b^2*(1-ws0) + b^2*d^2*ws0
        Risk <- list(asBias = b, asCov = asCov)
        a.fct <- function(x){ A * x[1] * z }
        body(a.fct) <- substitute({ A * x[1] * z }, list(A = A, z = z))
        Dom <- EuclideanSpace(dimension = 2)
        a <- EuclRandVarList(RealRandVariable(Map = list(a.fct), Domain = Dom))

        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info))    
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "Distribution",
                                          risk = "asBias", 
                                          neighbor = "Av2CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo, upper, maxiter, tol, warn){
        K <- E(Regressor, fun = function(x){ x %*% t(x) })
        z <- q(ErrorL2deriv)(0.5)
        Eu <- E(ErrorL2deriv, function(x, z){abs(x - z)}, z = z)
        b <- sqrt(sum(diag(trafo %*% solve(K) %*% t(trafo))))/Eu
        
        if(is(ErrorL2deriv, "AbscontDistribution")){
            ws0 <- 0
            d <- 0
        }else{
            ws0 <- d(ErrorL2deriv)(z)
            d <- (2*p(ErrorL2deriv)(z) - ws0 - 1)/ws0                
        }

        info <- c("minimum asymptotic bias (lower case) solution")
        asCov <- b^2*(1-ws0) + b^2*d^2*ws0
        Risk <- list(asBias = b, asCov = asCov)

        return(list(A = 1, z = z, b = b, d = d, risk = Risk, info = info))    
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "UnivariateDistribution",
                                          risk = "asBias", 
                                          neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo, upper, maxiter, tol, warn){
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        Ex <- E(Regressor, abs)
        b <- zi*as.vector(trafo)/(-m1df(ErrorL2deriv, 0)*Ex)
        p0 <- p(ErrorL2deriv)(0)
        if(is(ErrorL2deriv, "AbscontDistribution"))
            ws0 <- 0
        else
            ws0 <- d(ErrorL2deriv)(0)

        if(zi == 1){
            a0 <- -b*(1-p0)/(1-ws0)
            a.fct <- function(x){a0}
            body(a.fct) <- substitute({a0}, list(a0 = a0))
            a <- RealRandVariable(Map = list(a.fct), Domain = img(Regressor))
        }else{
            a0 <- b*(p0-ws0)/(1-ws0)
            a.fct <- function(x){a0}
            body(a.fct) <- substitute({a0}, list(a0 = a0))
            a <- RealRandVariable(Map = list(a.fct), Domain = img(Regressor))
        }
            
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asCov = a0^2*(p0-ws0) + (zi*a0+b)^2*(1-p0), asBias = b)

        return(list(A = A, a = a, b = b, d = 1, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "MultivariateDistribution",
                                          risk = "asBias", 
                                          neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, 
             Finfo, trafo, upper, maxiter, tol, warn){
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

            p <- nrow(trafo)
            k <- ncol(trafo)
            A <- matrix(erg$par[1:(p*k)], ncol=k, nrow=p)
            a <- erg$par[(p*k+1):length(erg$par)]

            # not yet implemented: computation of d
            d <- numeric(p)

            info <- c("minimum asymptotic bias (lower case) solution")
            Risk <- list(asBias = b)
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
            p <- nrow(trafo)
            A <- matrix(erg$par, ncol = ncol(trafo), nrow = p)
            a <- numeric(p)
            # not yet implemented: computation of d
            d <- numeric(p)

            info <- c("minimum asymptotic bias (lower case) solution")
            Risk <- list(asBias = b)        
        }
        
        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "MultivariateDistribution",
                                          risk = "asBias", 
                                          neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, 
             Finfo, trafo, upper, maxiter, tol, warn){
        abs.fct <- function(x, A){ 
            v <- t(x) %*% t(A)
            return(as.vector(sqrt(v %*% t(v))))
        }

        bmin.fct <- function(param, Regressor, trafo){
            A <- matrix(param, ncol = ncol(trafo), nrow = nrow(trafo))

            return(E(Regressor, abs.fct, A = A)/sum(diag(A %*% t(trafo))))
        }
        
        erg <- optim(as.vector(trafo), bmin.fct, method = "Nelder-Mead", 
                    control=list(reltol=tol, maxit=100*maxiter), Regressor = Regressor, trafo = trafo)
        z <- q(ErrorL2deriv)(0.5)
        b <- 1/(erg$value*E(ErrorL2deriv, function(x, z){abs(x - z)}, z = z))

        p <- nrow(trafo)
        A <- matrix(erg$par, ncol=ncol(trafo), nrow=p)

        k <- dimension(img(Regressor))
        a.fct <- function(x){ A %*% x[1:k]*z }
        body(a.fct) <- substitute({ A %*% x[1:k]*z }, list(k = k))
        Dom <- EuclideanSpace(dimension = k + 1)
        a <- EuclRandVarList(EuclRandVariable(Map = list(a.fct), Domain = Dom, 
                                              dimension = trunc(p)))
        
        # not yet implemented: computation of d
        d <- numeric(p)

        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asBias = b)

        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "MultivariateDistribution",
                                          risk = "asBias", 
                                          neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, RegSymm, 
             Finfo, trafo, upper, maxiter, tol, warn){
        abs.fct <- function(x, A){ 
            return(as.vector(sqrt(sum((A %*% x)^2))))
        }

        bmin.fct <- function(param, Regressor, trafo, Eu){
            A <- matrix(param, ncol = ncol(trafo), nrow = nrow(trafo))
            return(E(Regressor, abs.fct, A = A)/sum(diag(A %*% t(trafo))))
        }
        
        erg <- optim(as.vector(trafo), bmin.fct, method = "Nelder-Mead", 
                    control=list(reltol=tol, maxit=100*maxiter), Regressor = Regressor, trafo = trafo)
        b <- 1/(erg$value*(-m1df(ErrorL2deriv, 0)))

        p <- nrow(trafo)
        A <- matrix(erg$par, ncol=ncol(trafo), nrow=p)

        p0 <- p(ErrorL2deriv)(0)
        if(is(ErrorL2deriv, "AbscontDistribution"))
            ws0 <- 0
        else
            ws0 <- d(ErrorL2deriv)(0)

        a0 <- -b*(1-p0)/(1-ws0)
        a.fct <- function(x){a0}
        body(a.fct) <- substitute({a0}, list(a0 = a0))
        a <- RealRandVariable(Map = list(a.fct), Domain = img(Regressor))
            
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asCov = a0^2*(p0-ws0) + (a0+b)^2*(1-p0), asBias = b)
        
        return(list(A = A, a = a, b = b, d = 1, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "RealRandVariable",
                                          Regressor = "Distribution",
                                          risk = "asBias", 
                                          neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorSymm, RegSymm, 
             ErrorDistr, ErrorL2derivSymm, ErrorL2derivDistrSymm, 
             Finfo, trafo, upper, z.start, A.start, maxiter, tol, warn){
        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        abs.fctu <- function(u, xx, ErrorL2deriv, A, a0, k){ 
            L1 <- as.vector(evalRandVar(ErrorL2deriv, u))
            L2 <- c(xx*L1[1], L1[2:length(L1)])
            Y <- as.vector(A %*% L2 - a0)

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
        param <- erg$par
        p <- nrow(trafo)
        k <- ncol(trafo)
        A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
        a <- param[(p*k+1):length(param)]
        d <- numeric(p)
        # computation of 'd' not yet implemented!
        
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asBias = b)
        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "RealRandVariable",
                                          Regressor = "Distribution",
                                          risk = "asBias", 
                                          neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorSymm, RegSymm, 
             ErrorDistr, ErrorL2derivSymm, ErrorL2derivDistrSymm, 
             Finfo, trafo, upper, z.start, A.start, maxiter, tol, warn){
        if(is.null(z.start)) z.start <- function(x){numeric(1)}
        if(is.null(A.start)) A.start <- trafo
        
        stop("not yet implemented")

        return(list(A = A, a = a, b = b, d = 0*a, risk = Risk, info = info))
    })
