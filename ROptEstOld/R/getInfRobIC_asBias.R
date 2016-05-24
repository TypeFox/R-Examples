###############################################################################
## get minimum bias solutions
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asBias", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn){
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        z <- q(L2deriv)(0.5)
        b <- zi*as.vector(trafo)/E(L2deriv, function(x, z){abs(x - z)}, z = z)

        if(is(L2deriv, "AbscontDistribution"))
            ws0 <- 0
        else
            ws0 <- d(L2deriv)(z)
        if(ws0 > 0)
            d <- (2*p(L2deriv)(z) - ws0 - 1)/ws0
        else 
            d <- 0
        
        info <- c("minimum asymptotic bias (lower case) solution")
        asCov <- b^2*(1-ws0) + b^2*d^2*ws0
        Risk <- list(asBias = b, asCov = asCov)

        return(list(A = A, a = zi*z, b = b, d = d, risk = Risk, info = info))    
    })
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asBias", 
                                   neighbor = "TotalVarNeighborhood"),
    function(L2deriv, risk, neighbor, symm, Finfo, trafo, 
             upper, maxiter, tol, warn){
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        b <- zi*as.vector(trafo)/(-m1df(L2deriv, 0))
        p0 <- p(L2deriv)(0)
        if(is(L2deriv, "AbscontDistribution"))
            ws0 <- 0
        else
            ws0 <- d(L2deriv)(0)

        if(zi == 1)
            a <- -b*(1-p0)/(1-ws0)
        else
            a <- b*(p0-ws0)/(1-ws0)
            
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asCov = a^2*(p0-ws0) + (zi*a+b)^2*(1-p0), asBias = b)

        return(list(A = A, a = a, b = b, d = 1, risk = Risk, info = info))
    })
setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asBias", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, Finfo, z.start, A.start, trafo, upper, 
             maxiter, tol, warn){                
        if(is.null(z.start)) z.start <- numeric(ncol(trafo))
        if(is.null(A.start)) A.start <- trafo

        abs.fct <- function(x, L2, stand, cent){ 
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- apply(X, 2, "%*%", t(stand)) 

            return(sqrt(colSums(Y^2)))
        }
        bmin.fct <- function(param, L2deriv, Distr, trafo, z.comp){
            p <- nrow(trafo)
            k <- ncol(trafo)
            A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
            z <- numeric(k)
            z[z.comp] <- param[(p*k+1):length(param)]

            return(E(object = Distr, fun = abs.fct, L2 = L2deriv, stand = A, 
                     cent = z, useApply = FALSE)/sum(diag(A %*% t(trafo))))
        }
        
        nrvalues <- length(L2deriv)
        z.comp <- rep(TRUE, nrvalues)
        for(i in 1:nrvalues)
            if(is(L2derivDistrSymm[[i]], "SphericalSymmetry"))
                if(L2derivDistrSymm[[i]]@SymmCenter == 0)
                    z.comp[i] <- FALSE

        A.vec <- as.vector(A.start)
        erg <- optim(c(A.vec, z.start[z.comp]), bmin.fct, method = "Nelder-Mead", 
                    control = list(reltol = tol, maxit = 100*maxiter), 
                    L2deriv = L2deriv, Distr = Distr, trafo = trafo, z.comp = z.comp)
        b <- 1/erg$value
        param <- erg$par
        p <- nrow(trafo)
        k <- ncol(trafo)
        A <- matrix(param[1:(p*k)], ncol=k, nrow=p)
        z <- numeric(k)
        z[z.comp] <- param[(p*k+1):length(param)]
        a <- as.vector(A %*% z)
        d <- numeric(p)
        # computation of 'd', in case 'L2derivDistr' not abs. cont.
        
        info <- c("minimum asymptotic bias (lower case) solution")
        Risk <- list(asBias = b)

        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info))
    })
