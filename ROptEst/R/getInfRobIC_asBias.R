###############################################################################
## get minimum bias solutions
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asBias", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, symm, trafo, maxiter, 
             tol, warn, Finfo,
             verbose = NULL, ...){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
        erg <- minmaxBias(L2deriv = L2deriv, neighbor = neighbor,
                   biastype = biastype(risk), symm = symm, 
                   trafo = trafo, maxiter = maxiter, 
                   tol = tol, warn = warn, Finfo = Finfo)
        asCov <- erg$risk$asCov
        b <- erg$risk$asBias$value
        r <- neighbor@radius^2
        erg$risk <- c(erg$risk, 
                    list(trAsCov = list(value = asCov, normtype = NormType()),
                         asMSE = list(value = asCov + r^2*b^2, 
                                      r = r,
                                      at = neighbor)))
       return(erg)
    })
setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asBias", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, DistrSymm, L2derivSymm, 
             L2derivDistrSymm, z.start, 
             A.start, Finfo, trafo, maxiter, tol, warn,
             verbose = NULL, ...){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        k <- ncol(trafo); p <- nrow(trafo)
        if(is(neighbor,"TotalVarNeighborhood") && p>1)
           stop("Not yet implemented.")

        normtype <- normtype(risk)

#        if(FALSE){
        if(is(normtype,"SelfNorm")){
                warntxt <- paste(gettext(
                "Using self-standardization, there are problems with the existence\n"
                               ),gettext(
                "of a minmax Bias IC. Instead we compute the optimal MSE-solution\n"
                               ),gettext(
                "to a large radius (r = 15)\n"
                               ))
                if(warn) cat(warntxt)
                neighbor@radius <- 15
                res <- getInfRobIC(L2deriv = L2deriv, 
                        risk = asMSE(normtype = normtype), 
                        neighbor = neighbor, Distr = Distr, 
                        DistrSymm = DistrSymm, L2derivSymm = L2derivSymm, 
                        L2derivDistrSymm = L2derivDistrSymm, Finfo = Finfo, 
                        trafo = trafo, onesetLM = FALSE, z.start = z.start, 
                        A.start = A.start, upper = 1e4, maxiter = maxiter, 
                        tol = tol, warn = warn, verbose = verbose)
                
                A.max <- max(abs(res$A))
                res$A <- res$A/A.max
                res$a <- res$a/A.max
                w <- res$w
                A.max.w <- max(abs(stand(w)))
                stand(w) <- stand(w)/A.max.w
                normtype  <- res$normtype
                weight(w) <- minbiasweight(w, neighbor = neighbor,
                                           biastype = biastype(risk),
                                           normW = normtype)
                res$w <- w
                
                res$risk$asBias <- list(value = sqrt(nrow(trafo)), 
                                       biastype = symmetricBias(), 
                                       normtype = normtype, 
                                       neighbortype = class(neighbor),
                                       remark = gettext("value is only a bound"))
                return(res)
        }
#        }

        FI <- solve(trafo%*%solve(Finfo)%*%t(trafo))
        if(is(normtype,"QFNorm")) 
           {QuadForm(normtype) <- PosSemDefSymmMatrix(FI); 
            normtype(risk) <- normtype}

        ## determine which entries must be computed
        # by default everything
        z.comp <- rep(TRUE,k)
        A.comp <- matrix(rep(TRUE,k*k),nrow=k)

        # otherwise if trafo == unitMatrix may use symmetry info
        if(.isUnitMatrix(trafo)){
            comp <- .getComp(L2deriv, DistrSymm, L2derivSymm, L2derivDistrSymm)
            z.comp <- comp$"z.comp"
            A.comp <- comp$"A.comp"
        }

        return(minmaxBias(L2deriv = L2deriv, neighbor = neighbor,
                   biastype = biastype(risk), normtype = normtype(risk),
             Distr = Distr, z.start = z.start, A.start = A.start, 
             z.comp = z.comp, A.comp = A.comp, Finfo = Finfo, trafo = trafo,
             maxiter = maxiter, tol = tol, verbose = verbose))
    })



###############################################################################
## helper function minmaxBias
###############################################################################

setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "ContNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, symm, 
             trafo, maxiter, tol, warn, Finfo){
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
        Risk <- list(asBias = list(value = b, biastype = biastype, 
                                       normtype = NormType(), 
                                       neighbortype = class(neighbor)), 
                     asCov = asCov)

        w <- new("HampelWeight")
        cent(w) <- z
        stand(w) <- A
        clip(w) <- b
        weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biastype,
                       normW = NormType())
        return(list(A = A, a = zi*z, b = b, d = d, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType(),
                    problem = FALSE))
    })

setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "TotalVarNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, symm, trafo, 
             maxiter, tol, warn, Finfo){
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
            a <- -b*(p0-ws0)/(1-ws0)

        info <- c("minimum asymptotic bias (lower case) solution")
        asCov <- b^2/(1-ws0)*(1-p0)*(p0-ws0) #a^2*(p0-ws0) + (zi*a+b)^2*(1-p0)
        Risk <- list(asBias = list(value = b, biastype = biastype, 
                                       normtype = NormType(), 
                                       neighbortype = class(neighbor)), 
                     asCov = asCov)

        w <- new("BdStWeight")
        stand(w) <- A
        clip(w) <- c(a, a+b)
        weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biastype)
        return(list(A = A, a = a, b = b, d = 0, risk = Risk, info = info,
                    w = w, biastype = biastype, normtype = NormType(),
                    problem = FALSE))
    })

setMethod("minmaxBias", signature(L2deriv = "RealRandVariable", 
                                   neighbor = "ContNeighborhood", 
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, normtype, Distr, 
             z.start, A.start,  z.comp, A.comp, Finfo, trafo, maxiter,  tol,
             verbose = NULL){

        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
        DA.comp <- abs(trafo) %*% A.comp != 0
        eerg <- .LowerCaseMultivariate(L2deriv = L2deriv, neighbor = neighbor,
             biastype = biastype, normtype = normtype, Distr = Distr,
             Finfo = Finfo, trafo, z.start, A.start = A.start, z.comp = z.comp,
             A.comp = DA.comp, maxiter = maxiter, tol = tol, verbose = verbose)
        erg <- eerg$erg

        b <- 1/erg$value
        param <- erg$par
        lA.comp <- sum(DA.comp)
        
        p <- nrow(trafo)
        k <- ncol(trafo)
        A <- matrix(0, ncol=k, nrow=p)
        A[DA.comp] <- matrix(param[1:lA.comp], ncol=k, nrow=p)
        A.max <- max(abs(A))
        A <- A/A.max
        z <- numeric(k)
        z[z.comp] <- param[(lA.comp+1):length(param)]
        a <- as.vector(A %*% z)
        d <- numeric(p)

        # to be done: 
        # computation of 'd', in case 'L2derivDistr' not abs. cont.
        
        w <- eerg$w
        normtype <- eerg$normtype
        problem <- eerg$problem

        if(verbose)
           .checkPIC(L2deriv, neighbor, Distr, trafo, z, A, w, z.comp, A.comp)


        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype, Distr = Distr, 
                       V.comp = A.comp, cent = a, 
                       stand = A, w = w)

        std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)

        info <- c("minimum asymptotic bias (lower case) solution")
        trAsCov <- sum(diag(std%*%Cov))
        r <- neighbor@radius
        asMSE <- r^2 * b^2 + trAsCov
        Risk <- list(asBias = list(value = b, biastype = biastype, 
                                   normtype = normtype, 
                                   neighbortype = class(neighbor)), 
                     asCov = Cov,
                     trAsCov = list(value = trAsCov, normtype = normtype),
                     asMSE = list(value = asMSE, 
                                  r = r,
                                  at = neighbor))
        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = normtype,
                    problem = problem))
    })


setMethod("minmaxBias", signature(L2deriv = "RealRandVariable",
                                   neighbor = "TotalVarNeighborhood",
                                   biastype = "BiasType"),
    function(L2deriv, neighbor, biastype, normtype, Distr,
             z.start, A.start,  z.comp, A.comp, Finfo, trafo, maxiter,  tol,
             verbose = NULL){
        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        eerg <- .LowerCaseMultivariateTV(L2deriv = L2deriv,
             neighbor = neighbor, biastype = biastype,
             normtype = normtype, Distr = Distr, Finfo = Finfo, trafo = trafo,
             A.start = A.start, maxiter = maxiter,
             tol = tol, verbose = verbose)


        p <- nrow(trafo)
        k <- ncol(trafo)

        A <- eerg$A
        b <- eerg$b
        w <- eerg$w
        a <- eerg$a
        z <- numeric(k)
        d <- 0

        # to be done:
        # computation of 'd', in case 'L2derivDistr' not abs. cont.

        if(verbose)
           .checkPIC(L2deriv, neighbor, Distr, trafo, z, A, w,
                     z.comp=rep(TRUE,k), A.comp=matrix(TRUE,k,k))


        Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor,
                       biastype = biastype, Distr = Distr,
                       V.comp = matrix(TRUE), cent = numeric(k),
                       stand = A, w = w)

        std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)

        info <- c("minimum asymptotic bias (lower case) solution")
        trAsCov <- sum(diag(std%*%Cov))
        r <- neighbor@radius
        asMSE <- r^2 * b^2 + trAsCov
        Risk <- list(asBias = list(value = b, biastype = biastype,
                                   normtype = normtype,
                                   neighbortype = class(neighbor)),
                     asCov = Cov,
                     trAsCov = list(value = trAsCov, normtype = normtype),
                     asMSE = list(value = asMSE,
                                  r = r,
                                  at = neighbor))
        return(list(A = A, a = a, b = b, d = d, risk = Risk, info = info,
                    w = w, biastype = biastype, normtype = normtype,
                    problem = problem))
    })


setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "ContNeighborhood", 
                                   biastype = "asymmetricBias"),
    function(L2deriv, neighbor, biastype, symm, 
             trafo, maxiter, tol, warn, Finfo){
        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]
        zi <- sign(as.vector(trafo))
        A <- as.matrix(zi)
        z <- q(L2deriv)(nu1/(nu1+nu2))
        b <- zi*as.vector(trafo)/E(L2deriv, function(x, z){(x - z)*(x>z)/nu2 +
                 (z-x)*(z>x)/nu1}, z = z)

        b1 <- b / nu1
        b2 <- b / nu2

        p <- p(L2deriv)(z)

        if(is(L2deriv, "AbscontDistribution"))
            ws0 <- 0
        else
            ws0 <- d(L2deriv)(z)
        if(ws0 > 0)
            d <- (-b2*(1-p)+b1*(p-ws0))/ws0/b
        else
            d <- 0

        info <- c("minimum asymptotic bias (lower case) solution")
        asCov <- b2^2*(1-p)+b1^2*(p-ws0) + b^2*d^2*ws0
        Risk <- list(asBias = list(value = b, biastype = biastype, 
                                       normtype = NormType(), 
                                       neighbortype = class(neighbor)), 
                     asCov = asCov)

        w <- new("HampelWeight")
        cent(w) <- z
        stand(w) <- A
        clip(w) <- c(b1,b2)
        weight(w) <- minbiasweight(w, neighbor = neighbor, biastype = biastype)

        return(list(A = A, a = zi*z, b = b, d = d, risk = Risk, info = info, 
                    w = w, biastype = biastype, normtype = NormType(),
                    problem = FALSE))
           })

setMethod("minmaxBias", signature(L2deriv = "UnivariateDistribution", 
                                   neighbor = "ContNeighborhood", 
                                   biastype = "onesidedBias"),
    function(L2deriv, neighbor, biastype, symm, 
             trafo, maxiter, tol, warn, Finfo){
        infotxt <- c("minimum asymptotic bias (lower case) solution")
        noIC <- function(){
                warntxt <- paste(gettext(
                "There exists no IC attaining the infimal maxBias.\n"),
                                 gettext(
                "Instead we issue an IC with a very small Bias bound (starting with\n"), 
                                 gettext(
                "'tol'+ w_inf, w_inf = -1/inf_P psi or 1/sup_P psi).\n"
                                         ))
                w <- if(sign(biastype)>0) -1/q(L2deriv)(0) else 1/q(L2deriv)(1)
                if(warn) cat(warntxt)
                bd <- tol + w
                while (!is.list(try(
                res <- getInfRobIC(L2deriv = L2deriv, 
                        risk = asHampel(bound = bd, biastype = biastype), 
                        neighbor = neighbor, Finfo = Finfo, trafo = trafo, tol = tol,
                        warn = warn, noLow = TRUE, symm = symm, maxiter = maxiter,
                        verbose = verbose),
                silent = TRUE))) bd <- bd * 1.5
                return(res)}
        if(is(L2deriv, "DiscreteDistribution"))
           { if(is.finite(lowerCaseRadius(L2deriv, neighbor, risk = asMSE(), biastype)))
                {
                 sign <- sign(biastype)
                 w0 <- options("warn")
                 on.exit(options(w0))
                 options(warn = -1)
        
                 l <- length(support(L2deriv))
                 if (sign>0)
                      {z0 <- support(L2deriv)[1] 
                       deltahat <- support(L2deriv)[2]-z0
                 }else{
                       z0 <- support(L2deriv)[l]
                       deltahat <- z0-support(L2deriv)[l-1]
                 }
                 p0 <- d(L2deriv)(z0)   
                 v1 <- (1-p0)/p0/z0
                 v2 <- -1/z0
                 c0 <- deltahat*p0/2
                 A0 <- abs(1/z0/c0)
                 zc <- z0+sign(biastype)*deltahat*(1-p0)/2
                 a0 <- A0*zc
                 b0 <- abs(1/z0)
                 d0  <- 0 
                 asCov <- v1^2*(1-p0)+v2^2*p0
                 Risk0 <- list(asBias = list(value = b0, biastype = biastype, 
                               normtype = NormType(), 
                               neighbortype = class(neighbor)), 
                               asCov = asCov)
                 A0 <- matrix(A0,1,1)

                 w <- new("HampelWeight")
                 cent(w) <- z0
                 stand(w) <- A0
                 clip(w) <- b0
                 weight(w) <- minbiasweight(w, neighbor = neighbor, 
                               biastype = biastype)

                }else{return(noIC())}
            }else{return(noIC())}                    
        return(list(A = A0, a = a0, b = b0, d = d0, risk = Risk0, 
                    info = infotxt, w = w, biastype = biastype, 
                    normtype = NormType(),
                    problem = FALSE))
           })
