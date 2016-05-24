###############################################################################
## get optimally robust IC for convex asymptotic risks
###############################################################################
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "UnivariateDistribution", 
                                          risk = "asUnOvShoot", 
                                          neighbor = "UncondNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo, upper, maxiter, tol, warn){
        radius <- neighbor@radius
        if(identical(all.equal(radius, 0), TRUE)){
            if(warn) cat("'radius == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                        risk = asCov(), neighbor = TotalVarNeighborhood(), 
                        ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                        RegSymm = RegSymm, Finfo = Finfo, trafo = trafo)
            Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                   Regressor = Regressor, neighbor = neighbor,
                                   clip = res$b, cent = res$a, stand = res$A)
            res$risk <- c(Risk, res$risk)
            return(res)
        }

        bound <- risk@width*(-m1df(ErrorL2deriv, 0))*E(Regressor, abs)
        if(is(neighbor, "ContNeighborhood")){
            if(identical(all.equal(radius, 2*bound), TRUE)){
                zi <- sign(as.vector(trafo))
                A <- as.matrix(zi)
                b <- zi*as.vector(trafo)*2*risk@width/radius
                if(is(ErrorL2deriv, "AbscontDistribution"))
                    ws0 <- 0
                else
                    ws0 <- d(ErrorL2deriv)(0)
                if(is(Regressor, "AbscontDistribution"))
                    ws0x <- 0
                else
                    ws0x <- d(Regressor)(0)
                p0 <- ((1-p(Regressor)(0))*(p(ErrorL2deriv)(0)-ws0) 
                       + (p(Regressor)(0)-ws0x)*(1-p(ErrorL2deriv)(0)))
                p1 <- ((1-p(Regressor)(0))*(1-p(ErrorL2deriv)(0)) 
                       + (p(Regressor)(0)-ws0x)*(p(ErrorL2deriv)(0)-ws0))
                ws00 <- 1 - p0 - p1

                if(zi == 1)
                    a <- -b*p1/(1-ws00)
                else
                    a <- b*p0/(1-ws00)
            
                info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
                Risk <- list(asUnOvShoot = 0.5)

                return(list(A = A, a = a, b = b, d = 1, risk = Risk, info = info))                
            }
            if(radius > 2*bound)
                stop("boundedness condition is violated!")
        }

        if(is(neighbor, "TotalVarNeighborhood")){
            if(identical(all.equal(radius, bound), TRUE)){
                zi <- sign(as.vector(trafo))
                A <- as.matrix(zi)
                b <- zi*as.vector(trafo)*risk@width/radius
                if(is(ErrorL2deriv, "AbscontDistribution"))
                    ws0 <- 0
                else
                    ws0 <- d(ErrorL2deriv)(0)
                if(is(Regressor, "AbscontDistribution"))
                    ws0x <- 0
                else
                    ws0x <- d(Regressor)(0)
                p0 <- ((1-p(Regressor)(0))*(p(ErrorL2deriv)(0)-ws0) 
                       + (p(Regressor)(0)-ws0x)*(1-p(ErrorL2deriv)(0)))
                p1 <- ((1-p(Regressor)(0))*(1-p(ErrorL2deriv)(0)) 
                       + (p(Regressor)(0)-ws0x)*(p(ErrorL2deriv)(0)-ws0))
                ws00 <- 1 - p0 - p1

                if(zi == 1)
                    a <- -b*p1/(1-ws00)
                else
                    a <- b*p0/(1-ws00)
            
                info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
                Risk <- list(asUnOvShoot = 0.5)

                return(list(A = A, a = a, b = b, d = 1, risk = Risk, info = info))                
            }
            if(radius > bound)
                stop("boundedness condition is violated!")
        }

        z <- 0
        b <- 0
        if(is(ErrorL2derivDistrSymm, "SphericalSymmetry")) 
            z.comp <- !(ErrorL2derivDistrSymm@SymmCenter == 0)
        else
            z.comp <- TRUE

        iter <- 0
        repeat{
            iter <- iter + 1
            b.old <- b
            z.old <- z
            b <- try(uniroot(getInfClipRegTS, lower = .Machine$double.eps^0.75, 
                        upper = upper, tol = tol, ErrorL2deriv = ErrorL2deriv, 
                        Regressor = Regressor, risk = risk, neighbor = neighbor, 
                        z.comp = z.comp, cent = z)$root, silent = TRUE)

            if(!is.numeric(b)){
                if(warn) cat("Could not determine optimal clipping bound!\n", 
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
            z <- getInfCentRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor, 
                        neighbor = TotalVarNeighborhood(radius = neighbor@radius), 
                        clip = b, cent = z, z.comp = z.comp)

            prec <- max(abs(b-b.old), abs(z-z.old))
#            cat("current precision in IC algo:\t", prec, "\n")
            if(is(Regressor, "UnivariateDistribution") & (!z.comp)) break
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                    neighbor = TotalVarNeighborhood(), clip = b, cent = z)

        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                               Regressor = Regressor, neighbor = neighbor,
                               clip = b, cent = z, stand = A)

        return(list(A = as.matrix(A), a = A*z, b = A*b, d = NULL, risk = Risk, info = info))    
    })

setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution",
                                          Regressor = "UnivariateDistribution", 
                                          risk = "asUnOvShoot", 
                                          neighbor = "CondNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo, upper, maxiter, tol, warn){
        radiusCurve <- neighbor@radiusCurve
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

        radCx <- radiusCurve(x.vec)
        if(identical(all.equal(max(radCx), 0), TRUE)){
            if(warn) cat("'radiusCurve == 0' => (classical) optimal IC\n", 
                         "in sense of Cramer-Rao bound is returned\n")
            res <- getInfRobRegTypeIC(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                        risk = asCov(), neighbor = CondTotalVarNeighborhood(), 
                        ErrorL2derivDistrSymm = ErrorL2derivDistrSymm,
                        RegSymm = RegSymm, Finfo = Finfo, trafo = trafo)
            Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                                   Regressor = Regressor, neighbor = neighbor,
                                   clip = res$b, cent = res$a, stand = res$A)
            res$risk <- c(Risk, res$risk)
            return(res)
        }

        bound <- 1/(risk@width*(-m1df(ErrorL2deriv, 0)))
        if(is(neighbor, "CondContNeighborhood")){
            test <- abs(x.vec) - radCx*bound/2
            if(identical(all.equal(max(test), 0), TRUE) & any(test == 0)){
                if(!is(RegDistr, "AbscontDistribution"))
                    if(!identical(all.equal(d(Regressor)(0), 0), TRUE))
                        stop("Solution only available under 'K(x=0)!=0'!")

                stop("boundary case not yet implemented")
            }
            if(all(test < 0))
                stop("boundedness condition is violated!")
        }

        if(is(neighbor, "CondTotalVarNeighborhood")){
            test <- abs(x.vec) - radCx*bound
            if(identical(all.equal(max(test), 0), TRUE) & any(test == 0)){
                if(!is(RegDistr, "AbscontDistribution"))
                    if(!identical(all.equal(d(Regressor)(0), 0), TRUE))
                        stop("Solution only available under 'K(x=0)!=0'!")

                stop("boundary case not yet implemented")
            }
            if(all(test < 0))
                stop("boundedness condition is violated!")
        }

        z.vec <- numeric(length(x.vec))
        b.vec <- numeric(length(x.vec))
        if(is(ErrorL2derivDistrSymm, "SphericalSymmetry")) 
            z.comp <- !(ErrorL2derivDistrSymm@SymmCenter == 0)
        else
            z.comp <- TRUE

        iter <- 0
        repeat{
            iter <- iter + 1
            b.old <- b.vec
            z.old <- z.vec
            for(i in 1:length(x.vec)){
                if(test[i] <= 0){ # covers x == 0
                    b.vec[i] <- 0
                    z.vec[i] <- 0
                }else{ # here: x != 0
                    b.vec[i] <- try(uniroot(getInfClipRegTS, lower = .Machine$double.eps^0.75, 
                                    upper = upper, tol = tol, ErrorL2deriv = ErrorL2deriv, 
                                    Regressor = x.vec[i], risk = risk, neighbor = neighbor)$root, 
                                silent = TRUE)
                    z.vec[i] <- getInfCentRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = x.vec[i], 
                                neighbor = CondTotalVarNeighborhood(radius = neighbor@radius, 
                                                                    radiusCurve = neighbor@radiusCurve), 
                                clip = b.vec[i], cent = z.vec[i], z.comp = z.comp)
                }
            }

            prec <- max(abs(b.vec-b.old), abs(z.vec-z.old))
#            cat("current precision in IC algo:\t", prec, "\n")
            if(is(Regressor, "UnivariateDistribution") & (!z.comp)) break
            if(prec < tol) break
            if(iter > maxiter){
                cat("maximum iterations reached!\n", "achieved precision:\t", prec, "\n")
                break
            }
        }
        if(is(Regressor, "DiscreteDistribution")){
            bfun <- function(x){ 
                ind <- (round(x, 8) == round(x.vec, 8))
                if(any(ind))
                    return(b.vec[ind])
                else
                    return(NA)
            }
            zfun <- function(x){ 
                ind <- (round(x, 8) == round(x.vec, 8))
                if(any(ind))
                    return(z.vec[ind])
                else
                    return(NA)
            }
        }else{
            if(is.finite(q(Regressor)(0))){
                yleft.b <- NA
                yleft.z <- NA
            }else{
                yleft.b <- b.vec[1]
                yleft.z <- z.vec[1]
            }
            if(is.finite(q(Regressor)(1))){
                yright.b <- NA
                yright.z <- NA
            }else{
                yright.b <- b.vec[length(b.vec)]
                yright.z <- z.vec[length(z.vec)]
            }
            zfun <- approxfun(x.vec, z.vec, yleft = yleft.b, yright = yright.b)
            bfun <- approxfun(x.vec, b.vec, yleft = yleft.z, yright = yright.z)
        }

        A <- getInfStandRegTS(ErrorL2deriv = ErrorL2deriv, Regressor = Regressor,
                    neighbor = CondTotalVarNeighborhood(), clip = bfun, cent = zfun)

        info <- paste("optimally robust IC for", sQuote(class(risk)[1]))
        Risk <- getAsRiskRegTS(risk = risk, ErrorL2deriv = ErrorL2deriv, 
                               Regressor = Regressor, neighbor = neighbor,
                               clip = bfun, cent = zfun, stand = A)

        bfct <- function(x){bf <- bfun; A * bf(x)}
        body(bfct) <- substitute({bf <- bfun; A * bf(x)}, list(bfun = bfun, A = A))
        b <- RealRandVariable(Map = list(bfct), Domain = Reals())

        afct <- function(x){z <- zfun; A * z(x)}
        body(afct) <- substitute({z <- zfun; A * z(x)}, list(zfun = zfun, A = A))
        a <- RealRandVariable(Map = list(afct), Domain = Reals())

        return(list(A = as.matrix(A), a = a, b = b, d = NULL, risk = Risk, info = info))    
    })
