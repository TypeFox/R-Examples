###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "Distribution", 
                                          risk = "asGRisk", 
                                          neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo, upper, maxiter, tol, warn){
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                        risk = asCov(), neighbor = neighbor, 
                        ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                        RegSymm = RegSymm, Finfo = Finfo, trafo = trafo)
            Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                   Regressor = Regressor, neighbor = neighbor,
                                   clip = res$b, cent = res$a, stand = res$A, 
                                   trafo = trafo)
            res$risk <- c(Risk, res$risk)
            return(res)
        }
        z <- function(x){numeric(1)}

        A <- trafo
        b <- 0
        if(is(ErrorL2derivDistrSymm, "SphericalSymmetry")) 
            z.comp <- !(ErrorL2derivDistrSymm@SymmCenter == 0)
        else
            z.comp <- TRUE
        
        if(is(Regressor, "UnivariateDistribution")){
            if(is(Regressor, "AbscontDistribution")){
                xlower <- ifelse(is.finite(q(Regressor)(0)), q(Regressor)(0), q(Regressor)(distr::TruncQuantile))
                xupper <- ifelse(is.finite(q(Regressor)(1)), q(Regressor)(1), q(Regressor)(1 - distr::TruncQuantile))
                x.vec <- seq(from = xlower, to = xupper, length = 1000)
            }else{
                if(is(Regressor, "DiscreteDistribution"))
                    x.vec <- support(Regressor) 
                else
                    x.vec <- unique(r(Regressor)(distr::RtoDPQ.e))
            }
            z.vec <- numeric(length(x.vec))
        }else{
            if(is(Regressor, "DiscreteMVDistribution"))
                x.vec <- support(Regressor)
            else
                stop("not yet implemented")
    
            z.vec <- numeric(nrow(x.vec))
        }

        iter <- 0
        repeat{
            iter <- iter + 1
            b.old <- b
            z.vec.old <- z.vec
            A.old <- A
            b <- try(uniroot(getInfClipRegTS, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, ErrorL2deriv = ErrorL2deriv, 
                        Regressor = Regressor, risk = risk, neighbor = neighbor, 
                        z.comp = z.comp, stand = A, cent = z)$root, silent = TRUE)
            if(!is.numeric(b)){
                if(warn) cat("Could not determine optimal clipping bound!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n")
                res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                                risk = asBias(), neighbor = neighbor, 
                                ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
                                trafo = trafo, maxiter = maxiter, tol = tol, warn = warn)
                Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                       Regressor = Regressor, neighbor = neighbor,
                                       clip = res$b, cent = res$a, stand = res$A, 
                                       trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            zz <- getInfCentRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                        neighbor = neighbor, clip = b, cent = z, stand = A, 
                        z.comp = z.comp, x.vec = x.vec)
            z.vec <- zz$z.vec
            z <- zz$z.fct
            if(is(Regressor, "UnivariateDistribution"))
                A <- A.old
            else
                A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                            neighbor = neighbor, z.comp = z.comp, clip = b, cent = z, 
                            stand = A, trafo = trafo)

            prec <- max(abs(b-b.old), abs(A-A.old), abs(z.vec-z.vec.old))
            if(is(Regressor, "MultivariateDistribution") & z.comp)
                cat("current precision in IC algo:\t", prec, "\n")
            if(is(Regressor, "UnivariateDistribution") & (!z.comp)) break
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }

        if(is(Regressor, "UnivariateDistribution")){
            A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                        neighbor = neighbor, z.comp = z.comp, clip = b, cent = z, 
                        stand = A, trafo = trafo)
            b <- as.vector(A) * b
        }
        p <- dimension(img(Regressor))
        a.fct <- function(x){ A %*% x[1:p]*z(x[1:p]) }
        body(a.fct) <- substitute({ A %*% x[1:p]*z(x[1:p]) }, list(p = p))
        Dom <- EuclideanSpace(dimension = p + 1)
        a <- EuclRandVarList(EuclRandVariable(Map = list(a.fct), Domain = Dom, 
                                              dimension = trunc(nrow(trafo))))

        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                               Regressor = Regressor, neighbor = neighbor,
                               clip = b, cent = a, stand = A, trafo = trafo)
        Risk <- c(Risk, list(asBias = b))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "RealRandVariable",
                                          Regressor = "Distribution", 
                                          risk = "asGRisk", 
                                          neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorSymm, RegSymm, 
             ErrorDistr, ErrorL2derivSymm, ErrorL2derivDistrSymm, 
             Finfo, trafo, upper, z.start, A.start, maxiter, tol, warn){
        k <- length(ErrorL2deriv)
        if(is.null(z.start)) z.start <- function(x){numeric(k)}
        if(is.null(A.start)) A.start <- trafo

        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                        risk = asCov(), neighbor = neighbor, ErrorDistr = ErrorDistr, 
                        Finfo = Finfo, trafo = trafo)
            Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                   Regressor = Regressor, neighbor = neighbor,
                                   clip = res$b, cent = res$a, stand = res$A, 
                                   trafo = trafo)
            res$risk <- c(Risk, res$risk)
            return(res)
        }

        z <- z.start
        if(is(Regressor, "UnivariateDistribution")){
            if(is(Regressor, "AbscontDistribution")){
                xlower <- ifelse(is.finite(q(Regressor)(0)), q(Regressor)(0), q(Regressor)(distr::TruncQuantile))
                xupper <- ifelse(is.finite(q(Regressor)(1)), q(Regressor)(1), q(Regressor)(1 - distr::TruncQuantile))
                x.vec <- seq(from = xlower, to = xupper, length = 1000)
            }else{
                if(is(Regressor, "DiscreteDistribution"))
                    x.vec <- support(Regressor) 
                else
                    x.vec <- unique(r(Regressor)(distr::RtoDPQ.e))
            }
            z.vec <- matrix(0, ncol = k, nrow = length(x.vec))
        }else{
            if(is(Regressor, "DiscreteMVDistribution"))
                x.vec <- support(Regressor)
            else
                stop("not yet implemented")

            z.vec <- matrix(0, ncol = k, nrow = nrow(x.vec))
        }

        A <- A.start
        b <- 0
        nrvalues <- length(ErrorL2deriv)
        z.comp <- rep(TRUE, nrvalues)
        A.comp <- matrix(TRUE, ncol = nrvalues, nrow = nrvalues)
        for(i in 1:nrvalues){
            if(is(ErrorL2derivDistrSymm[[i]], "SphericalSymmetry"))
                if(ErrorL2derivDistrSymm[[i]]@SymmCenter == 0)
                    z.comp[i] <- FALSE
        }
        for(i in 1:(nrvalues-1))
            for(j in (i+1):nrvalues){
                if(is(ErrorSymm, "SphericalSymmetry")){
                    if((is(ErrorL2derivSymm[[i]], "OddSymmetric") & is(ErrorL2derivSymm[[j]], "EvenSymmetric"))
                       | (is(ErrorL2derivSymm[[j]], "OddSymmetric") & is(ErrorL2derivSymm[[i]], "EvenSymmetric")))
                        if((ErrorL2derivSymm[[i]]@SymmCenter == ErrorL2derivSymm[[j]]@SymmCenter)
                           & (ErrorL2derivSymm[[i]]@SymmCenter == ErrorSymm@SymmCenter))
                            A.comp[i,j] <- FALSE
                }
            }
        A.comp[col(A.comp) < row(A.comp)] <- A.comp[col(A.comp) > row(A.comp)]

        iter <- 0
        repeat{
            iter <- iter + 1
            z.old <- z
            z.vec.old <- z.vec
            b.old <- b
            A.old <- A
            tb <- system.time(b <- try(uniroot(getInfClipRegTS, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, ErrorL2deriv = ErrorL2deriv, 
                        Regressor = Regressor, risk = risk, neighbor = neighbor, 
                        ErrorDistr = ErrorDistr, stand = A, cent = z, 
                        trafo = trafo)$root, silent = FALSE))
            cat("time for b-step:\t", tb, "\n")
            if(!is.numeric(b)){
                if(warn) cat("Could not determine optimal clipping bound!\n", 
                             "'radius >= maximum radius' for the given risk?\n",
                             "=> the minimum asymptotic bias (lower case) solution is returned\n",
                             "If 'no' => Try again with modified starting values ",
                             "'z.start' and 'A.start'\n") 
                res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                                risk = asBias(), neighbor = neighbor, ErrorSymm = ErrorSymm,
                                RegSymm = RegSymm, ErrorDistr = ErrorDistr, 
                                ErrorL2derivSymm = ErrorL2derivSymm, 
                                ErrorL2derivDistrSymm = ErrorL2derivDistrSymm, 
                                Finfo = Finfo, trafo = trafo, z.start = z.start, A.start = A.start,
                                upper = upper, maxiter = maxiter, tol = tol, warn = warn)
                Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                       Regressor = Regressor, neighbor = neighbor,
                                       clip = res$b, cent = res$a, stand = res$A, 
                                       trafo = trafo)
                res$risk <- c(Risk, res$risk)
                return(res)
            }
            tz <- system.time(zz <- getInfCentRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                            neighbor = neighbor, ErrorDistr = ErrorDistr, stand = A, 
                            cent = z, clip = b, z.comp = z.comp, x.vec = x.vec))
            cat("time for z-step:\t", tz, "\n")
            z.vec <- zz$z.vec
            z <- zz$z.fct

            tA <- system.time(A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                        neighbor = neighbor, ErrorDistr = ErrorDistr, A.comp = A.comp,
                        stand = A, clip = b, cent = z, trafo = trafo))
            cat("time for A-step:\t", tA, "\n")
            prec <- max(abs(b-b.old), max(abs(A-A.old)), max(abs(z.vec-z.vec.old)))
            cat("current precision in IC algo:\t", prec, "\n")
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        p <- dimension(img(Regressor))
        a.fct <- function(x){
            z <- z(x[1:p])
            A %*% c(x[1:p]*z[1], z[2:k])
        }
        body(a.fct) <- substitute({z <- z(x[1:p])
                                   A %*% c(x[1:p]*z[1], z[2:k])},
                                  list(k = length(ErrorL2deriv), p = p))
        Dom <- EuclideanSpace(dimension = p + 1)
        a <- EuclRandVarList(EuclRandVariable(Map = list(a.fct), Domain = Dom, dimension = p + k - 1))

        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                               Regressor = Regressor, neighbor = neighbor,
                               clip = b, cent = a, stand = A, trafo = trafo)
        Risk <- c(Risk, list(asBias = b))

        return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))
    })
