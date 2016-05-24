###############################################################################
## radius minimax optimally robust IC 
## for L2ParamFamily and asymptotic risks
###############################################################################
setMethod("radiusMinimaxIC", signature(L2Fam = "L2ParamFamily", 
                                       neighbor = "UncondNeighborhood",
                                       risk = "asGRisk"),
    function(L2Fam, neighbor, risk, loRad = 0, upRad = Inf, z.start = NULL,
             A.start = NULL, upper = NULL, lower = NULL,
             OptOrIter = "iterate", maxiter = 50,
             tol = .Machine$double.eps^0.4, warn = FALSE,
             verbose = NULL, loRad0 = 1e-3, ..., returnNAifProblem = FALSE,
             loRad.s = NULL, upRad.s = NULL){
        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
        ow <- options("warn")
        on.exit(options(ow))
        if(length(loRad) != 1)
            stop("'loRad' is not of length == 1")
        if(length(upRad) != 1)
            stop("'upRad' is not of length == 1")
        if(loRad >= upRad)
            stop("'upRad < loRad' is not fulfilled")

        biastype <- biastype(risk)
        L2derivDim <- numberOfMaps(L2Fam@L2deriv)
        trafo <- trafo(L2Fam@param)

        if(is(normtype(risk),"SelfNorm")||is(normtype(risk),"InfoNorm"))
           upRad <- min(upRad,10) 

        ineff <- NULL
        args.IC <- list(L2deriv = NULL, neighbor = NULL,
                       risk = risk, symm = NULL,
                       Finfo = L2Fam@FisherInfo, upper = upper, lower = lower,
                       trafo = trafo, maxiter = maxiter, tol = tol,
                       warn = warn, verbose = verbose)

        args.R  <- list(risk = risk, L2deriv = NULL,
                         neighbor = neighbor, biastype = biastype,
                         normtype = normtype(risk), trafo = trafo)

        args.Ie <- list(radius = NULL, L2Fam = L2Fam, neighbor,
                        risk = risk, upper.b = upper, lower.b = lower,
                        loRad = loRad, upRad = upRad, eps = tol,
                        MaxIter = maxiter, warn = warn,
                        loNorm = NormType(), upNorm = NormType(),
                        verbose=verbose, withRetIneff = TRUE)
        fct.Ie <- function(x){
                 args.Ie$radius  <- x
#                 print(with(args.Ie, list(loRisk,upRisk,loRad,upRad)))
                 res <- do.call(getIneffDiff,args.Ie)
#                 print(res)
                 ineff <<- res["up"]
                 res["up"] - res["lo"]
        }

        .getRisk <- function(rad, fac = 1){
               neighbor@radius <- rad
               args.IC$maxiter <- maxiter * fac
               args.IC$neighbor <- args.R$neighbor <- neighbor
               res <- do.call(getInfRobIC, args.IC)
               args.R$clip <- res$b; args.R$stand <- res$A; args.R$cent <- res$a
               Norm <- NormType()
               if(L2derivDim > 1){
                   risk.0 <- risk; normtype(risk.0) <- Norm <- res$normtype
                   args.R$risk <- risk.0
               }
               return(list(Risk=do.call(getAsRisk, args.R)[[1]], Norm=Norm))
        }

        if(L2derivDim == 1){
            options(warn = -1)
            args.R$L2deriv <- args.IC$L2deriv <- L2Fam@L2derivDistr[[1]]
            args.IC$symm <-  L2Fam@L2derivDistrSymm[[1]]

            if(is(neighbor,"TotalVarNeighborhood")) upper <- upper/2

            args.Ie$loRisk <- if(identical(all.equal(loRad, 0), TRUE))
                1/as.vector(L2Fam@FisherInfo) else .getRisk(loRad, 6)$Risk

            if(upRad == Inf){
                args.lR <- args.R
                args.lR$risk <- asBias(biastype = biastype)
                args.Ie$upRisk <- (do.call(getAsRisk, args.lR)$asBias)^2
            }else args.Ie$upRisk <- .getRisk(upRad)$Risk
#            print(c(rlo=loRad, Rlo=args.Ie$loRisk, rup=upRad,Rup=args.Ie$upRisk))
        }else{
            if(is(L2Fam@distribution, "UnivariateDistribution")){
               L2derivSymm <- L2Fam@L2derivSymm
               L2derivDistrSymm <- L2Fam@L2derivDistrSymm
               if((length(L2Fam@L2deriv) == 1) &
                     is(L2Fam@L2deriv[[1]], "RealRandVariable")){
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
               li.0 <- list(z.start = z.start, A.start = A.start,
                            OptOrIter = OptOrIter)
               li.1 <- list(Distr = L2Fam@distribution,
                            DistrSymm = L2Fam@distrSymm,
                            L2derivSymm = L2derivSymm,
                            L2derivDistrSymm = L2derivDistrSymm)

               args.R$L2deriv <- args.IC$L2deriv <- L2deriv
               args.IC <- c(args.IC, li.1, li.0)
               args.Ie <- c(args.Ie, li.0)

               normtype <- normtype(risk)
               Finfo <- L2Fam@FisherInfo

               p <- nrow(trafo)
               FI0 <- trafo%*%solve(Finfo)%*%t(trafo)

               if(is(normtype,"InfoNorm") || is(normtype,"SelfNorm") )
                    {QuadForm(normtype) <- PosSemDefSymmMatrix(solve(FI0));
                     normtype(risk) <- normtype}
               std <- if(is(normtype,"QFNorm")) QuadForm(normtype) else diag(p)
               loRisk <- sum(diag(std%*%FI0))

               options(warn = -1)

               if(identical(all.equal(loRad, 0), TRUE)){
                   loRad <- 0
                   args.Ie$loRisk <- loRisk
                   args.Ie$loNorm <- normtype
               }else{
                   rL <- .getRisk(loRad)
                   args.Ie$loRisk <- rL$Risk; args.Ie$loNorm <- rL$Norm
               }
               if(upRad == Inf){
                   args.lR <- c(list(risk = asBias(biastype = biastype(risk),
                                                 normtype = normtype),
                                L2deriv = L2deriv, neighbor = neighbor,
                                biastype = biastype, normtype = normtype(risk)),
                                li.1, list(Finfo = Finfo, trafo = trafo,
                                z.start = z.start, A.start = A.start,
                                maxiter = maxiter, tol = tol,
                                warn = warn, verbose = verbose))
                   biasR <- do.call(getAsRisk, args.lR)
                   args.Ie$upNorm <- biasR$normtype
                   args.Ie$upRisk <- (biasR$asBias)^2
               }else{
                   rL <- .getRisk(upRad)
                   args.Ie$upRisk <- rL$Risk; args.Ie$upNorm <- rL$Norm
               }

            }else{
                stop("not yet implemented")
            }
        }
        
        lower <- if(is.null(loRad.s)) max(loRad, loRad0) else loRad.s
        upper <- if(is.null(upRad.s)) {
             if(upRad == Inf) max(lower+2, 4) else upRad } else upRad.s
        leastFavR <- try(
                    uniroot(fct.Ie, lower = lower, upper = upper,
                         tol = .Machine$double.eps^0.25)$root , silent = TRUE)

        if(is(leastFavR, "try-error")){
           if(returnNAifProblem) return(NA)
           warnRund <- 1; isE <- TRUE
           fl <- (0.2/lower)^(1/6); fu <- (0.5/upper)^(1/6)
           while(warnRund < 7 && isE ){
              warnRund <- warnRund + 1
              lower <- lower * fl;  upper <- upper *fu
              if(is.finite(upRad)){
                 args.Ie$upRad <- upper; rL <- .getRisk(upper)
                 args.Ie$upRisk <- rL$Risk; args.Ie$upNorm <- rL$Norm
              }
              if(loRad>0){
                 args.Ie$loRad <- lower; rL <- .getRisk(lower)
                 args.Ie$loRisk <- rL$Risk; args.Ie$loNorm <- rL$Norm
              }
              leastFavR <- try(
                         uniroot(fct.Ie, lower = lower, upper = upper,
                         tol = .Machine$double.eps^0.25)$root, silent = TRUE)
              isE <- is(leastFavR, "try-error")
              if(isE) print(conditionMessage(attr(leastFavR,"condition")))
           }
           if(isE) stop("Problem: Zero search in getIneffDiff did not converge.")
           else warning(paste("Had to modify radius bounds to [", lower,
                        upper, "] after", warnRund, "iterations."))
        }
        neighbor@radius <- leastFavR
        args.IC$neighbor <- args.R$neighbor <- neighbor
        args.IC$returnNAifProblem <- returnNAifProblem
        res <- do.call(getInfRobIC, args.IC)
        if(returnNAifProblem) if(!is.null(res$problem)) if(res$problem) return(NA)
        options(ow)
        res$info <- c("radiusMinimaxIC", paste("radius minimax IC for radius interval [",
                        round(loRad, 3), ", ", round(upRad, 3), "]", sep=""))
        res$info <- rbind(res$info, c("radiusMinimaxIC",
                        paste("least favorable radius: ", round(leastFavR, 3), sep="")))
        res$info <- rbind(res$info, c("radiusMinimaxIC",
                        paste("maximum ", sQuote(class(risk)[1]), "-inefficiency: ",
                        round(ineff, 3), sep="")))
        res <- c(res, modifyIC = getModifyIC(L2FamIC = L2Fam,
                                             neighbor = neighbor,
                                             risk = risk))
        return(generateIC(neighbor, L2Fam, res))
            })
