###############################################################################
## asymptotic MSE
###############################################################################           
setMethod("getAsRisk", signature(risk = "asMSE",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "Neighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL, stand, trafo, ...){
        if(!is.finite(neighbor@radius))
            mse <- Inf
        else
            mse <- as.vector(stand)*as.vector(trafo)
        return(list(asMSE = mse))
    })

setMethod("getAsRisk", signature(risk = "asMSE",
                                 L2deriv = "EuclRandVariable",
                                 neighbor = "Neighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL, stand, trafo, ...){
        if(!is.finite(neighbor@radius))
            mse <- Inf
        else{
            p <- nrow(trafo)
            std <- if(is(normtype(risk),"QFNorm")) 
                      QuadForm(normtype(risk)) else diag(p)
                        
            mse <- sum(diag( std %*% stand %*% t(trafo)))
        }
        return(list(asMSE = mse))
    })

###############################################################################
## minimum asymptotic Bias
###############################################################################
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL, stand = NULL, trafo, ...){
        z <- q(L2deriv)(0.5)                                
        bias <- abs(as.vector(trafo))/E(L2deriv, function(x, z){abs(x - z)}, 
                                        useApply = FALSE, z = z)

        return(list(asBias = bias))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "TotalVarNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL, stand = NULL, trafo, ...){
        bias <- abs(as.vector(trafo))/(-m1df(L2deriv, 0))

        return(list(asBias = bias))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, 
             clip = NULL, cent = NULL, stand = NULL, Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, Finfo, trafo, z.start, A.start,  maxiter, tol, warn,
             verbose = NULL, ...){
        
        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        normtype <- normtype(risk)
        biastype <- biastype(risk)


        if(is(normtype,"SelfNorm")){
                warntxt <- paste(gettext(
                "Using self-standardization, there are problems with the existence\n"
                               ),gettext(
                "of a minmax Bias IC. Instead we return a lower bound.\n"
                               ))
                if(warn) cat(warntxt)
                return(list(asBias = sqrt(nrow(trafo))))        
        }
        comp <- .getComp(L2deriv, DistrSymm, L2derivSymm, L2derivDistrSymm)
        z.comp <- comp$"z.comp"
        A.comp <- comp$"A.comp"
        DA.comp <- abs(trafo) %*% A.comp != 0
        
        eerg <- .LowerCaseMultivariate(L2deriv = L2deriv, neighbor = neighbor, 
             biastype = biastype, normtype = normtype, Distr = Distr,  Finfo = Finfo,
             trafo = trafo, z.start = z.start, A.start = A.start, z.comp = z.comp,
             A.comp = DA.comp,  maxiter = maxiter, tol = tol, verbose = verbose)
        erg <- eerg$erg
        bias <- 1/erg$value
        
        return(list(asBias = bias, normtype = eerg$normtype))
    })
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "TotalVarNeighborhood",
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL,
             clip = NULL, cent = NULL, stand = NULL,  Distr, DistrSymm, L2derivSymm,
             L2derivDistrSymm, Finfo, trafo, z.start, A.start,  maxiter, tol, warn,
             verbose = NULL, ...){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        normtype <- normtype(risk)
        biastype <- biastype(risk)
        Finfo <- FisherInfo(L2deriv)

        eerg <- .LowerCaseMultivariateTV(L2deriv = L2deriv,
             neighbor = neighbor, biastype = biastype,
             normtype = normtype, Distr = Distr, Finfo = Finfo, trafo = trafo,
             A.start = A.start, maxiter = maxiter,
             tol = tol, verbose = verbose)
        erg <- eerg$b

        return(list(asBias = bias, normtype = eerg$normtype))
    })


###############################################################################
## asymptotic covariance
###############################################################################
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL,clip, cent, stand, trafo = NULL, ...){
#        c0 <- clip/abs(as.vector(stand))
#        D1 <- L2deriv - cent/as.vector(stand)
#        Cov <- (clip^2*(p(D1)(-c0) + 1 - p(D1)(c0))
#                + as.vector(stand)^2*(m2df(D1, c0) - m2df(D1, -c0)))
        Cov <- getInfV(L2deriv, neighbor, biastype(risk), clip/abs(as.vector(stand)), 
                cent/abs(as.vector(stand)), abs(as.vector(stand)))

        if(!is.null(trafo)) Cov <- trafo%*%Cov%*%t(trafo)
        return(list(asCov = Cov ))
    })
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "TotalVarNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand, trafo = NULL, ...){
         if(missing(biastype)||is.null(biastype)) biastype <- biastype(risk)
#        g0 <- cent/abs(as.vector(stand))
#        c0 <- clip/abs(as.vector(stand))
#        Cov <- (abs(as.vector(stand))^2*(g0^2*p(L2deriv)(g0) 
#                + (g0+c0)^2*(1 - p(L2deriv)(g0+c0))
#                + m2df(L2deriv, g0+c0) - m2df(L2deriv, g0)))
#        return(list(asCov = Cov))
        Cov  <-         getInfV(L2deriv, neighbor, biastype, clip/abs(as.vector(stand)), 
                cent/abs(as.vector(stand)), abs(as.vector(stand)))
        if(!is.null(trafo)) Cov <- trafo%*%Cov%*%t(trafo)
        return(list(asCov = Cov))

    })
setMethod("getAsRisk", signature(risk = "asCov",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip=NULL, cent, 
             stand, Distr, trafo = NULL, V.comp =  matrix(TRUE, ncol = nrow(stand), nrow = nrow(stand)), 
             w, ...){
        if(missing(biastype)||is.null(biastype)) biastype <- biastype(risk)
        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype, Distr = Distr, 
                       V.comp = V.comp, cent = cent, 
                       stand = stand, w = w)
        if(!is.null(trafo)) Cov <- trafo%*%Cov%*%t(trafo)
        return(list(asCov = Cov))
        })
#        Y <- as(stand %*% L2deriv - cent, "EuclRandVariable")
#        absY <- sqrt(Y %*% Y)
#        
#        nrvalues <- nrow(stand)
#        ICfct <- vector(mode = "list", length = nrvalues)
#        for(i in 1:nrvalues){
#            ICfct[[i]] <- function(x){ Yi(x)*pmin(1, b/absY(x)) }
#            body(ICfct[[i]]) <- substitute({ Yi(x)*pmin(1, b/absY(x)) },
#                                    list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = clip))
#        }
#        IC <- RealRandVariable(Map = ICfct, Domain = Y@Domain, Range = Y@Range)
#        Cov <- matrix(E(Distr, IC %*% t(IC)), ncol = nrvalues)
#
#        return(list(asCov = Cov))
#    })




###############################################################################
## trace of asymptotic covariance
###############################################################################
setMethod("getAsRisk", signature(risk = "trAsCov",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand, trafo=NULL, ...){
        if(missing(biastype)||is.null(biastype)) biastype <- biastype(risk)
        if(missing(normtype)||is.null(normtype)) normtype <- normtype(risk)
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype, clip = clip, cent = cent, stand = stand)$asCov
        std <- if(is(normtype,"QFNorm"))
                  QuadForm(normtype) else diag(nrow(as.matrix(Cov)))

        return(list(trAsCov = sum(diag(std%*%as.matrix(Cov)))))
    })

setMethod("getAsRisk", signature(risk = "trAsCov",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype, clip, cent, stand, Distr, trafo = NULL,
             V.comp =  matrix(TRUE, ncol = nrow(stand), nrow = nrow(stand)), w,...){
        if(missing(biastype)||is.null(biastype)) biastype <- biastype(risk)
        if(missing(normtype)||is.null(normtype)) normtype <- normtype(risk)
        Cov <- getAsRisk(risk = asCov(), L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype, Distr = Distr, clip = clip,
                         cent = cent, stand = stand, trafo = trafo,
                         V.comp =  V.comp, w = w)$asCov

        p <- nrow(stand)
        std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)
        
        return(list(trAsCov = sum(diag(std%*%Cov))))
    })

###############################################################################
## Anscombe criterion
###############################################################################
setMethod("getAsRisk", signature(risk = "asAnscombe",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip, cent, stand, trafo = NULL, FI, ... ){
        if(missing(biastype)||is.null(biastype)) biastype <- biastype(risk)
        if(missing(normtype)||is.null(normtype)) normtype <- normtype(risk)
        trAsCov.0 <- getAsRisk(risk = trAsCov(), L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype, normtype = normtype,
                         clip = clip, cent = cent, stand = stand)$trAsCov
        return(list(asAnscombe = FI/trAsCov.0))
    })
setMethod("getAsRisk", signature(risk = "asAnscombe",
                                 L2deriv = "RealRandVariable",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype, clip, cent, stand, Distr, trafo = NULL, 
             V.comp =  matrix(TRUE, ncol = nrow(stand), nrow = nrow(stand)), FI, 
             w, ...){
        if(missing(biastype)||is.null(biastype)) biastype <- biastype(risk)
        if(missing(normtype)||is.null(normtype)) normtype <- normtype(risk)
        trAsCov.0 <-  getAsRisk(risk = trAsCov(), L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype, normtype = normtype, 
                         Distr = Distr, clip = clip,  
                         cent = cent, stand = stand, V.comp = V.comp, 
                         w = w)$trAsCov
        return(list(asAnscombe = FI/trAsCov.0))
    })

###############################################################################
## asymptotic under-/overshoot risk
###############################################################################
setMethod("getAsRisk", signature(risk = "asUnOvShoot",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "UncondNeighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL,clip, cent, stand, trafo, ...){
        if(identical(all.equal(neighbor@radius, 0), TRUE))
            return(list(asUnOvShoot = pnorm(-risk@width/sqrt(as.vector(stand)))))
        
        g0 <- cent/abs(as.vector(stand))
        c0 <- clip/abs(as.vector(stand))
        s <- sqrt(g0^2*p(L2deriv)(g0) 
                  + (g0+c0)^2*(1 - p(L2deriv)(g0+c0))
                  + m2df(L2deriv, g0+c0) - m2df(L2deriv, g0))

        return(list(asUnOvShoot = pnorm(-risk@width*s)))
    })

###############################################################################
## asymptotic onesided bias
###############################################################################
setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "onesidedBias"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL,
             clip = NULL, cent = NULL, stand = NULL, trafo, ...){

        D1 <- L2deriv
        if(!is(D1, "DiscreteDistribution")) 
            return(list(asBias = 0, warn = gettext("not attained by IC")))

        sign <- sign(biastype)
        w0 <- options("warn")
        on.exit(options(w0))
        options(warn = -1)
        
        l <- length(support(L2deriv))
        if (sign>0)
           {z0 <- support(L2deriv)[1]; deltahat <- support(L2deriv)[2]-z0}
        else
           {z0 <- support(L2deriv)[l]; deltahat <- z0-support(L2deriv)[l-1]}

        bias <- abs(as.vector(trafo))/abs(z0)
        return(list(asBias = bias))
    })

###############################################################################
## asymptotic asymmetric bias
###############################################################################

setMethod("getAsRisk", signature(risk = "asBias",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "ContNeighborhood", 
                                 biastype = "asymmetricBias"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL,
             clip = NULL, cent = NULL, stand = NULL, trafo, ...){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]
        num <- nu2/(nu1+nu2)        
        z <- q(L2deriv)(num)
        Int <- E(L2deriv, function(x, m){abs(x-m)}, m = z)
        omega <- 2/(Int/nu1+Int/nu2)
        bias <- abs(as.vector(trafo))*omega
        return(list(asBias = bias))
    })

###############################################################################
## asymptotic semivariance
###############################################################################

setMethod("getAsRisk", signature(risk = "asSemivar",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "Neighborhood", 
                                 biastype = "onesidedBias"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL,
             clip, cent, stand, trafo, ...){
        A <- as.vector(stand)*as.vector(trafo)
        r <- neighbor@radius
        b <- clip*A
        if (sign(biastype)>0)
            v <- E(L2deriv, function(x) A^2*pmin(x-cent,clip)^2)
        else
            v <- E(L2deriv, function(x) A^2*pmax(x-cent,-clip)^2)
        sv <- r*b/sqrt(v)
        if(!is.finite(r))
            semvar <- Inf
        else
            semvar <- (v+r^2*b^2)*pnorm(sv)+ r*b*sqrt(v)*dnorm(sv)
        return(list(asSemivar = semvar))
    })

###############################################################################
## asymptotic L1 risk
###############################################################################           
setMethod("getAsRisk", signature(risk = "asL1",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "Neighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL, stand, trafo, ...){
        if(!is.finite(neighbor@radius))
            L1 <- Inf
        else{
             s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand)
             r <- neighbor@radius
             L1 <- get.asGRisk.fct(risk)(r,s=s^.5,b=clip)
        }
        return(list(asL1 = L1))
    })

###############################################################################
## asymptotic L4 risk
###############################################################################           
setMethod("getAsRisk", signature(risk = "asL4",
                                 L2deriv = "UnivariateDistribution",
                                 neighbor = "Neighborhood", 
                                 biastype = "ANY"),
    function(risk, L2deriv, neighbor, biastype, normtype = NULL, clip = NULL, cent = NULL, stand, trafo, ...){
        if(!is.finite(neighbor@radius))
            L4 <- Inf
        else{
             s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand)
             r <- neighbor@radius
             L4 <- get.asGRisk.fct(risk)(r,s=s^.5,b=clip)
        }
        return(list(asL4 = L4))
    })
