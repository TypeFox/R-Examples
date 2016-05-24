###############################################################################
## asymptotic MSE
###############################################################################
setMethod("getAsRisk", signature(risk = "asMSE",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "Neighborhood"),
    function(risk, L2deriv, neighbor, clip, cent, stand, trafo){
        if(!is.finite(neighbor@radius))
            mse <- Inf
        else
            mse <- as.vector(stand)*as.vector(trafo)
        return(list(asMSE = mse))
    })
setMethod("getAsRisk", signature(risk = "asMSE",
                                 L2deriv = "EuclRandVariable",
                                 neighbor = "Neighborhood"),
    function(risk, L2deriv, neighbor, clip, cent, stand, trafo){
        if(!is.finite(neighbor@radius))
            mse <- Inf
        else
            mse <- sum(diag(stand %*% t(trafo)))
        return(list(asMSE = mse))
    })

###############################################################################
## minimum asymptotic Bias
###############################################################################
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood"),
    function(risk, L2deriv, neighbor, trafo){
        z <- q(L2deriv)(0.5)
        bias <- abs(as.vector(trafo))/E(L2deriv, function(x, z){abs(x - z)}, 
                                        useApply = FALSE, z = z)

        return(list(asBias = bias))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "TotalVarNeighborhood"),
    function(risk, L2deriv, neighbor, trafo){
        bias <- abs(as.vector(trafo))/(-m1df(L2deriv, 0))

        return(list(asBias = bias))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood"),
    function(risk, L2deriv, neighbor, Distr, L2derivDistrSymm, trafo, 
             z.start, A.start,  maxiter, tol){                
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
        bias <- 1/erg$value
        
        return(list(asBias = bias))
    })

###############################################################################
## asymptotic covariance
###############################################################################
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood"),
    function(risk, L2deriv, neighbor, clip, cent, stand){
        c0 <- clip/abs(as.vector(stand))
        D1 <- L2deriv - cent/as.vector(stand)
        Cov <- (clip^2*(p(D1)(-c0) + 1 - p(D1)(c0))
                + as.vector(stand)^2*(m2df(D1, c0) - m2df(D1, -c0)))

        return(list(asCov = Cov))
    })
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "TotalVarNeighborhood"),
    function(risk, L2deriv, neighbor, clip, cent, stand){
        g0 <- cent/abs(as.vector(stand))
        c0 <- clip/abs(as.vector(stand))
        Cov <- (abs(as.vector(stand))^2*(g0^2*p(L2deriv)(g0) 
                + (g0+c0)^2*(1 - p(L2deriv)(g0+c0))
                + m2df(L2deriv, g0+c0) - m2df(L2deriv, g0)))

        return(list(asCov = Cov))
    })
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood"),
    function(risk, L2deriv, neighbor, Distr, clip, cent, stand){
        Y <- as(stand %*% L2deriv - cent, "EuclRandVariable")
        absY <- sqrt(Y %*% Y)
        
        nrvalues <- nrow(stand)
        ICfct <- vector(mode = "list", length = nrvalues)
        for(i in 1:nrvalues){
            ICfct[[i]] <- function(x){ Yi(x)*pmin(1, b/absY(x)) }
            body(ICfct[[i]]) <- substitute({ Yi(x)*pmin(1, b/absY(x)) },
                                    list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = clip))
        }
        IC <- RealRandVariable(Map = ICfct, Domain = Y@Domain, Range = Y@Range)
        Cov <- matrix(E(Distr, IC %*% t(IC)), ncol = nrvalues)

        return(list(asCov = Cov))
    })

###############################################################################
## trace of asymptotic covariance
###############################################################################
setMethod("getAsRisk", signature(risk = "trAsCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood"),
    function(risk, L2deriv, neighbor, clip, cent, stand){
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
                         clip = clip, cent = cent, stand = stand)$asCov

        return(list(trAsCov = as.vector(Cov)))
    })
setMethod("getAsRisk", signature(risk = "trAsCov",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood"),
    function(risk, L2deriv, neighbor, Distr, clip, cent, stand){
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
                         Distr = Distr, clip = clip, cent = cent, stand = stand)$asCov

        return(list(trAsCov = sum(diag(Cov))))
    })

###############################################################################
## asymptotic under-/overshoot risk
###############################################################################
setMethod("getAsRisk", signature(risk = "asUnOvShoot",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood"),
    function(risk, L2deriv, neighbor, clip, cent, stand, trafo){
        if(identical(all.equal(neighbor@radius, 0), TRUE))
            return(list(asUnOvShoot = pnorm(-risk@width/sqrt(as.vector(stand)))))
        
        g0 <- cent/abs(as.vector(stand))
        c0 <- clip/abs(as.vector(stand))
        s <- sqrt(g0^2*p(L2deriv)(g0) 
                  + (g0+c0)^2*(1 - p(L2deriv)(g0+c0))
                  + m2df(L2deriv, g0+c0) - m2df(L2deriv, g0))

        return(list(asUnOvShoot = pnorm(-risk@width*s)))
    })
