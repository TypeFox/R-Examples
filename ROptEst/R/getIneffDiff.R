###############################################################################
## get inefficiency difference for asymptotic MSE
###############################################################################
setMethod("getIneffDiff", signature(radius = "numeric", 
                                    L2Fam = "L2ParamFamily", 
                                    neighbor = "UncondNeighborhood",
                                    risk = "asMSE"),
    function(radius, L2Fam, neighbor, risk, loRad, upRad, loRisk, upRisk, 
             z.start = NULL, A.start = NULL, upper.b = NULL, lower.b = NULL,
             OptOrIter = "iterate", MaxIter, eps, warn,
             loNorm = NULL, upNorm = NULL,
             verbose = NULL, ..., withRetIneff = FALSE){
        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")

        L2derivDim <- numberOfMaps(L2Fam@L2deriv)
        if(L2derivDim == 1){
            ##print(radius)
            neighbor@radius <- radius
            res <- getInfRobIC(L2deriv = L2Fam@L2derivDistr[[1]], neighbor = neighbor, 
                        risk = risk, symm = L2Fam@L2derivDistrSymm[[1]], 
                        Finfo = L2Fam@FisherInfo, upper = upper.b, lower = lower.b,
                        trafo = trafo(L2Fam@param),
                        maxiter = MaxIter, tol = eps,
                        warn = warn, verbose = verbose)
            trafo <- as.vector(trafo(L2Fam@param))
            ineffLo <- (as.vector(res$A)*trafo - res$b^2*(radius^2-loRad^2))/loRisk
            ####cat("---------------\n")
            ##res00=res;res00$w <- NULL; res00$biastype <- NULL; res00$d <- NULL
            ##res00$normtype <- NULL;res00$info <- NULL;res00$risk <- NULL;
            ##print(res00)
            ##print(c(lower.b,upper.b,loRisk,"upR"=upRisk))
            ####cat("---------------\n")
            if(upRad == Inf)
                ineffUp <- res$b^2/upRisk
            else
                ineffUp <- (as.vector(res$A)*trafo - res$b^2*(radius^2-upRad^2))/upRisk
            ##assign("ineff", ineffUp, envir = sys.frame(which = -5))
            ##print(c(ineffUp,ineffLo,ineffUp - ineffLo))
             if(withRetIneff) return(c(lo= ineffLo, up=ineffUp))
             else return(ineffUp - ineffLo)
        }else{
            if(is(L2Fam@distribution, "UnivariateDistribution")){
                if((length(L2Fam@L2deriv) == 1) & is(L2Fam@L2deriv[[1]], "RealRandVariable")){
                    L2deriv <- L2Fam@L2deriv[[1]]
                    L2derivSymm <- L2Fam@L2derivSymm
                    L2derivDistrSymm <- L2Fam@L2derivDistrSymm
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
                trafo <- trafo(L2Fam@param)
                p <- nrow(trafo)
                neighbor@radius <- radius
                res <- getInfRobIC(L2deriv = L2deriv, neighbor = neighbor, risk = risk, 
                            Distr = L2Fam@distribution, DistrSymm = L2Fam@distrSymm, 
                            L2derivSymm = L2derivSymm, L2derivDistrSymm = L2derivDistrSymm, 
                            Finfo = L2Fam@FisherInfo, trafo = trafo, z.start = z.start, 
                            A.start = A.start, upper = upper.b, lower = lower.b,
                            OptOrIter = OptOrIter, maxiter = MaxIter,
                            tol = eps, warn = warn, verbose = verbose,
                            withPICcheck = FALSE,...)
                normtype(risk) <- res$normtype
                std <- if(is(normtype(risk),"QFNorm"))
                          QuadForm(normtype(risk)) else diag(p)

                biasLo <- biasUp <- res$b

                if(is(normtype(risk),"SelfNorm")){
                   IC <- generateIC(neighbor = neighbor, L2Fam = L2Fam, res = res)
                   biasLoE <- getBiasIC(IC = as(IC, "IC"), neighbor = neighbor, L2Fam = L2Fam, 
                                       biastype = symmetricBias(),
                                       normtype = loNorm, tol = eps, 
                                       numbeval = 1e4)
                   biasLo <- biasLoE$asBias$value
                   biasUpE <- getBiasIC(IC = as(IC, "IC"), neighbor = neighbor, L2Fam = L2Fam, 
                                       biastype = symmetricBias(),
                                       normtype = upNorm, tol = eps, 
                                       numbeval = 1e4)
                   biasUp <- biasUpE$asBias$value
                   ineffLo <- (p+biasLo^2*loRad^2)/loRisk
                   ineffUp <- if(upRad == Inf) biasUp^2/upRisk else
                                   (p+biasUp^2*upRad^2)/upRisk
                }else{
                    ineffLo <- (sum(diag(std%*%res$A%*%t(trafo))) - 
                                biasLo^2*(radius^2-loRad^2))/loRisk
                   if(upRad == Inf)
                      ineffUp <- biasUp^2/upRisk
                   else
                      ineffUp <- (sum(diag(std%*%res$A%*%t(trafo))) -
                                  biasUp^2*(radius^2-upRad^2))/upRisk
                }
                if(verbose)
                    cat(paste(rep("-",75), sep = "", collapse = ""),"\n",
                        "current radius:   ", round(radius,4),
                        "\tMSE-inefficiency difference:   ",
                        round(ineffUp - ineffLo,4),
                        paste("\n",paste(rep("-",75), sep = "",
                                         collapse = ""),"\n",sep="")
                        )

                if(withRetIneff) return(c(lo= ineffLo, up=ineffUp))
                else return(ineffUp - ineffLo)
            }else{
                stop("not yet implemented")
            }
        }
    })
