###############################################################################
## radius minimax optimally robust IC
## for L2ParamFamily and asymptotic risks
###############################################################################
setMethod("leastFavorableRadius", signature(L2Fam = "L2ParamFamily",
                                            neighbor = "UncondNeighborhood",
                                            risk = "asGRisk"),
    function(L2Fam, neighbor, risk, rho, upRad = 1,
            z.start = NULL, A.start = NULL, upper = 100,
            OptOrIter = "iterate", maxiter = 100,
            tol = .Machine$double.eps^0.4, warn = FALSE, verbose = NULL){
        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
        if(length(rho) != 1)
            stop("'rho' is not of length == 1")
        if((rho <= 0)||(rho >= 1))
            stop("'rho' not in (0,1)")

        biastype <- biastype(risk)
        normtype <- normtype(risk)

        trafo <- trafo(L2Fam@param)
        FI0 <- trafo%*%solve(L2Fam@FisherInfo)%*%t(trafo)
        FI <- solve(FI0)
        if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") )
           {QuadForm(normtype) <- PosSemDefSymmMatrix(FI);
            normtype(risk) <- normtype}

        L2derivDim <- numberOfMaps(L2Fam@L2deriv)
        if(L2derivDim == 1){
            leastFavFct <- function(r, L2Fam, neighbor, risk, rho,
                                    upper.b, MaxIter, eps, warn){
                loRad <- r*rho
                upRad <- r/rho
                lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
                upper <- ifelse(upRad == Inf, 10, upRad)
                ow <- options("warn")
                on.exit(options(ow))
                options(warn = -1)
                .getRisk <- function(rad, fac = 1){
                    neighbor@radius <- rad
                    res <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]],
                                neighbor = neighbor,
                                risk = risk, symm = L2Fam@L2derivDistrSymm[[1]],
                                Finfo = L2Fam@FisherInfo, upper = upper.b,
                                trafo = trafo, maxiter = MaxIter*fac, tol = eps,
                                warn = warn, verbose = verbose)
                   return(getAsRisk(risk = risk, L2deriv = L2Fam@L2derivDistr[[1]],
                                        neighbor = neighbor, biastype = biastype,
                                        normtype = normtype,
                                        clip = res$b, cent = res$a,
                                        stand = res$A, trafo = trafo)[[1]])
                 }

                if(identical(all.equal(loRad, 0), TRUE)){
                    loRad <- 0
                    loRisk <- 1/as.vector(L2Fam@FisherInfo)
                }else{
                    loRisk <- .getRisk(loRad,6)
                }

                if(upRad == Inf){
                    bmin <- getAsRisk(risk = asBias(biastype = biastype),
                                L2deriv = L2Fam@L2derivDistr[[1]],
                                neighbor = neighbor, biastype = biastype,
                                normtype = normtype,
                                trafo = trafo, symm = L2Fam@L2derivSymm[[1]])
                    upRisk <- bmin^2
                }else{
                    upRisk <- .getRisk(upRad)
                }
                loNorm<- upNorm <- NormType()
                args.Ie <- list(L2Fam = L2Fam, neighbor = neighbor, risk = risk,
                                loRad = loRad, upRad = upRad, loRisk = loRisk,
                                upRisk = upRisk, upper.b = upper.b, eps = eps,
                                MaxIter = MaxIter, warn = warn,
                                loNorm = loNorm, upNorm = upNorm,
                                withRetIneff = TRUE)
                ineff <- NULL
                fct.Ie <- function(x){
                  args.Ie$radius  <- x
                  res <- do.call(getIneffDiff,args.Ie)
                  ineff <<- res["up"]
                  res["up"] - res["lo"]
                }
                leastFavR <- try(uniroot(fct.Ie, lower = lower, upper = upper,
                                 tol = .Machine$double.eps^0.25)$root, silent =TRUE)
                if(is(leastFavR, "try-error")){
                   warnRund <- 1; isE <- TRUE
                   while(warnRund < 7 && isE ){
                     warnRund <- warnRund + 1
                     lower <- lower * 2;  upper <- upper / 2
                     if( warnRund == 4 ) min(upper, 1.5)
                     if(is.finite(upRad)){
                        args.Ie$upRad <- upper; args.Ie$upRisk <- .getRisk(upper)
                     }
                     if(loRad>0){
                        args.Ie$loRad <- lower; args.Ie$loRisk <- .getRisk(lower)
                     }
                     leastFavR <- try(
                         uniroot(fct.Ie, lower = lower, upper = upper,
                         tol = .Machine$double.eps^0.25)$root, silent = TRUE)
                     isE <- is(leastFavR, "try-error")
                     if(isE) print(conditionMessage(attr(leastFavR,"condition")))
                   }
                   if(isE)
                      stop("Problem: Zero search in getIneffDiff did not converge.")
                   else warning(paste("Had to modify radius bounds to [", lower,
                        upper, "] after", warnRund, "iterations."))
                }

                options(ow)
                cat("current radius:\t", r, "\tinefficiency:\t", ineff, "\n")
                return(ineff)
            }
            leastFavR <- optimize(leastFavFct, lower = 1e-4, upper = upRad,
                            tol = .Machine$double.eps^0.25, maximum = TRUE,
                            L2Fam = L2Fam, neighbor = neighbor, risk = risk,
                            rho = rho, upper.b = upper, MaxIter = maxiter,
                            eps = tol, warn = warn)

            res <- list(rho = rho, leastFavorableRadius = leastFavR$maximum,
                        ineff = leastFavR$objective)
            names(res)[3] <- paste(class(risk)[1], "-inefficiency", sep="")
            return(res)
        }else{
            if(is(L2Fam@distribution, "UnivariateDistribution")){
                L2derivSymm <- L2Fam@L2derivSymm
                L2derivDistrSymm <- L2Fam@L2derivDistrSymm
                if((length(L2Fam@L2deriv) == 1) & is(L2Fam@L2deriv[[1]], "RealRandVariable")){
                    L2deriv <- L2Fam@L2deriv[[1]]
                }else{
                    L2deriv <- diag(dimension(L2Fam@L2deriv)) %*% L2Fam@L2deriv
                    L2deriv <- RealRandVariable(Map = L2deriv@Map, Domain = L2deriv@Domain)
                    nrvalues <- numberOfMaps(L2deriv)
                    if(numberOfMaps(L2Fam@L2deriv) != nrvalues){
                        L1 <- vector("list", nrvalues)
                        L2 <- vector("list", nrvalues)
                        for(i in 1:nrvalues){
                            L1[[i]] <- NonSymmetric()
                            L2[[i]] <- NoSymmetry()
                        }
                        L2derivSymm <- new("FunSymmList", L1)
                        L2derivDistrSymm <- new("DistrSymmList", L2)
                    }
                }

                std <- if(is(normtype,"QFNorm"))
                       QuadForm(normtype) else diag(nrow(trafo))


                leastFavFct <- function(r, L2Fam, neighbor, risk, rho,
                                        z.start, A.start, upper.b, MaxIter, eps, warn){

                    loRad <- r*rho
                    upRad <- r/rho
                    lower <- ifelse(identical(all.equal(loRad, 0), TRUE), 1e-4, loRad)
                    upper <- ifelse(upRad == Inf, 10, upRad)
                    ow <- options("warn")
                    on.exit(options(ow))
                    options(warn = -1)

                    .getRisk <- function(rad, fac = 1){
                        neighbor@radius <- rad
                        res <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor,
                                        risk = risk, Distr = L2Fam@distribution,
                                        DistrSymm = L2Fam@distrSymm,
                                        L2derivSymm = L2derivSymm,
                                        L2derivDistrSymm = L2derivDistrSymm,
                                        Finfo = L2Fam@FisherInfo, trafo = trafo,
                                        z.start = z.start, A.start = A.start,
                                        upper = upper.b, OptOrIter = OptOrIter,
                                        maxiter = MaxIter, tol = eps, warn = warn,
                                        verbose = verbose)
                        risk.0 <- risk
                        normtype(risk.0) <- .Norm <- res$normtype
                        .Risk <- getAsRisk(risk = risk.0, L2deriv = L2deriv,
                                          neighbor = neighbor, biastype = biastype,
                                          normtype = normtype, clip = res$b,
                                          cent = res$a, stand = res$A,
                                          trafo = trafo)[[1]]
                        return(list(Risk=.Risk,Norm=.Norm))
                    }

                    if(identical(all.equal(loRad, 0), TRUE)){
                        loRad <- 0
                        loRisk <- sum(diag(std%*%FI0))
                        loNorm <- normtype
                    }else{
                        rL <- .getRisk(loRad)
                        loRisk <- rL$Risk
                        loNorm <- rL$Norm
                    }

                    if(upRad == Inf){
                        biasR <- getAsRisk(risk = asBias(biastype = biastype(risk),
                                      normtype = normtype), L2deriv = L2deriv,
                                      neighbor = neighbor, biastype = biastype,
                                      normtype = normtype,
                                      Distr = L2Fam@distribution,
                                      DistrSymm = L2Fam@distrSymm,
                                      L2derivSymm = L2derivSymm,
                                      L2derivDistrSymm= L2derivDistrSymm,
                                Finfo = L2Fam@FisherInfo, trafo = trafo,
                                z.start = z.start, A.start = A.start,
                                maxiter = maxiter, tol = tol,
                                warn = warn, verbose = verbose)
                        bmin <- biasR$asBias
                        upRisk <- bmin^2
                        upNorm <- biasR$normtype
                    }else{
                        rL <- .getRisk(upRad)
                        upRisk <- rL$Risk
                        upNorm <- rL$Norm
                    }
                    args.Ie <- list(L2Fam = L2Fam,
                                    neighbor = neighbor, z.start = z.start,
                                    A.start = A.start, upper.b = upper.b,
                                    risk = risk,
                                    loRad = loRad, upRad = upRad,
                                    loRisk = loRisk, upRisk = upRisk,
                                    eps = eps, OptOrIter = OptOrIter,
                                    MaxIter = MaxIter, warn = warn,
                                    loNorm = loNorm, upNorm = upNorm,
                                    withRetIneff = TRUE)
                    ineff <- NULL
                    fct.Ie <- function(x){
                      args.Ie$radius  <- x
                      res <- do.call(getIneffDiff,args.Ie)
                      ineff <<- res["up"]
                      res["up"] - res["lo"]
                    }
                    leastFavR <- try(uniroot(fct.Ie, lower = lower, upper = upper,
                                 tol = .Machine$double.eps^0.25)$root, silent =TRUE)
                    if(is(leastFavR, "try-error")){
                       warnRund <- 1; isE <- TRUE
                       fl <- (0.2/lower)^(1/6); fu <- (0.5/upper)^(1/6)
                       while(warnRund < 7 && isE ){
                         warnRund <- warnRund + 1
                         lower <- lower * fl;  upper <- upper *fr
                         if( warnRund == 4 ) min(upper, 1.5)
                         if(is.finite(upRad)){
                            args.Ie$upRad <- upper; rL <- .getRisk(upper)
                            args.Ie$upRisk <- rL$Risk; args.Ie$upNorm <- rL$Norm
                         }
                         if(loRad>0){
                            args.Ie$upRad <- upper; rL <- .getRisk(upper)
                            args.Ie$upRisk <- rL$Risk; args.Ie$upNorm <- rL$Norm
                         }
                         leastFavR <- try(
                             uniroot(fct.Ie, lower = lower, upper = upper,
                             tol = .Machine$double.eps^0.25)$root, silent = TRUE)
                         isE <- is(leastFavR, "try-error")
                         if(isE) print(conditionMessage(attr(leastFavR,"condition")))
                       }
                       if(isE)
                          stop("Problem: Zero search in getIneffDiff did not converge.")
                       else warning(paste("Had to modify radius bounds to [", lower,
                            upper, "] after", warnRund, "iterations."))
                    }
                    options(ow)

                    if(verbose)
                       cat(paste(rep("-",75), sep = "", collapse = ""),"\n")
                    cat("current radius:   ", round(r,4),
                        "\tinefficiency:   ", round(ineff,4))
                    if(verbose)
                       cat(paste("\n",paste(rep("-",75), sep = "",
                                        collapse = ""),"\n",sep=""))
                    else cat("\n")

                    return(ineff)
                }
                if(is.null(z.start)) z.start <- numeric(L2derivDim)
                if(is.null(A.start)) A.start <- trafo
                leastFavR <- optimize(leastFavFct, lower = 1e-4, upper = upRad,
                                tol = .Machine$double.eps^0.25, maximum = TRUE,
                                L2Fam = L2Fam, neighbor = neighbor, risk = risk,
                                rho = rho, z.start = z.start, A.start = A.start,
                                upper.b = upper, MaxIter = maxiter, eps = tol, warn = warn)

                res <- list(rho = rho, leastFavorableRadius = leastFavR$maximum,
                            ineff = leastFavR$objective)
                names(res)[3] <- paste(class(risk)[1], "-inefficiency", sep="")
                return(res)
            }else{
                stop("not yet implemented")
            }
        }
    })
