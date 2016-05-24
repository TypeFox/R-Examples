###############################################################################
## one-step estimator
###############################################################################

oneStepEstimator <- function(x, IC, start = NULL,
                             useLast = getRobAStBaseOption("kStepUseLast"),
                             withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                             IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                             na.rm = TRUE, startArgList = NULL, ...){
        es.call <- match.call()
        es.call[[1]] <- as.name("oneStepEstimator")

        if(! is(IC, "InfluenceCurve"))
           stop("Argument 'IC' must be of class 'InfluenceCurve'")

        if(is(IC, "IC")){
            erg <- kStepEstimator(x = x, IC = IC, start = start, steps = 1L,
                           useLast = useLast, withUpdateInKer = withUpdateInKer,
                           IC.UpdateInKer = IC.UpdateInKer, na.rm = na.rm,
                           startArgList = startArgList, ...)
            Infos(erg) <- gsub("kStep","oneStep", Infos(erg))
            erg@estimate.call <- es.call
            return(erg)
        }

        ### now: IC is not of class "IC" -> have to do it by hand:

        if(withUpdateInKer || !is.null(IC.UpdateInKer))
           warning("We do not use args 'withUpdateInKer' or 'IC.UpdateInKer' in case arg 'IC' is not of class 'IC'")

        ### transform if necessary
        x0 <- x
        x0 <- if(is.numeric(x) && ! is.matrix(x)) {
                 x0 <- as.matrix(x)
              }
        completecases <- complete.cases(x0)
        if(na.rm) x0 <- na.omit(x0)

        if(missing(start)||is.null(start))
           stop("In case arg 'IC' is not of class 'IC', arg 'start' must not be missing.")

        nrvalues <- dimension(IC@Curve)
        start.val <- kStepEstimator.start(start, x=x0, nrvalues = nrvalues, na.rm = na.rm,
                                          L2Fam = NULL, startList = startArgList)

        res <- start.val + rowMeans(evalIC(IC, x0), na.rm = na.rm)
        nms <- if(!is.null(dim(res))) colnames(res)  else names(start.val)
        dim(res) <- NULL
        if(is.null(names(res))) names(res) <- nms
        Infos <- matrix(c("oneStepEstimator", "1-step estimate"), ncol = 2)
        colnames(Infos) <- c("method", "message")
        asVar <- NULL
        asBias <- NULL

        nuis.idx <- if(is(start,"Estimate")) start@nuis.idx else NULL
        fixed <- if(is(start,"Estimate")) start@fixed else NULL

        new("kStepEstimate", name = "1-step estimate", estimate = res,
            untransformed.estimate = res, untransformed.asvar = NULL,
             fixed = fixed, nuis.idx = nuis.idx,
             completecases = completecases,
            estimate.call = es.call, samplesize = nrow(x0), asvar = asVar,
            asbias = asBias, pIC = IC, steps = 1L, Infos = Infos,
            start = start, startval = start.val, ustartval = start.val)
    }




### old Routine
#setMethod("oneStepEstimator", signature(x = "numeric",
#                                        IC = "InfluenceCurve",
#                                        start = "ANY"),
#    function(x, IC, start, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("oneStepEstimator")
#        nrvalues <- dimension(IC@Curve)
#        if(is(start, "Estimate")){
#            start0 <- estimate(start)
#        }else{
#            start0 <- start
#        }
#        if(is.list(start0)) start0 <- unlist(start0)
#        if(nrvalues != length(start0))
#            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
#
#        res <- start0 + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)
#
#        if(is(IC, "IC")){
#            L2Fam <- eval(CallL2Fam(IC))
#            Infos <- matrix(c("oneStepEstimator",
#                            paste("1-step estimate for", name(L2Fam))),
#                            ncol = 2)
#            colnames(Infos) <- c("method", "message")
#            if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("oneStepEstimator",
#                                        "computation of IC, asVar and asBias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("oneStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("oneStepEstimator",
#                                        "computation of IC, asVar and asBias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#        }else{
#            Infos <- matrix(c("oneStepEstimator", "1-step estimate"), ncol = 2)
#            colnames(Infos) <- c("method", "message")
#            asVar <- NULL
#            asBias <- NULL
#        }
#        trafo <- trafo(L2Fam, newParam)$mat
#        if(.isUnitMatrix(trafo)){
#           uest <- res
#           uvar <- asVar
#        }else{
#           st0 <- if("untransformed.estimate" %in% getSlots(class(start)))
#                  start@untransformed.estimate
#                  else ginv(trafo)%*%start0
#           uest <- st0 + ginv(trafo)%*%(res-start0)
#           uvar <- NULL
#        }
#
#        new("kStepEstimate", name = "1-step estimate", estimate = res,
#            estimate.call = es.call, samplesize = length(x), asvar = asVar,
#            untransformed.estimate = uest, untransformed.asvar = uvar,
#            asbias = asBias, pIC = IC, steps = 1L, Infos = Infos)
#    })
#setMethod("oneStepEstimator", signature(x = "matrix",
#                                        IC = "InfluenceCurve",
#                                        start = "ANY"),
#    function(x, IC, start, useLast = getRobAStBaseOption("kStepUseLast")){
#        es.call <- match.call()
#        es.call[[1]] <- as.name("oneStepEstimator")
#        nrvalues <- dimension(IC@Curve)
#        if(is(start, "Estimate")){
#            start0 <- estimate(start)
#        }else{
#            start0 <- start
#        }
#        if(is.list(start0)) start0 <- unlist(start0)
#        if(nrvalues != length(start0))
#            stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
#        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
#            stop("'x' has wrong dimension")
#
#        res <- start0 + rowMeans(evalIC(IC, x), na.rm = TRUE)
#
#        if(is(IC, "IC")){
#            L2Fam <- eval(CallL2Fam(IC))
#            Infos <- matrix(c("oneStepEstimator",
#                            paste("1-step estimate for", name(L2Fam))),
#                            ncol = 2)
#            colnames(Infos) <- c("method", "message")
#            if(is(L2Fam, "L2GroupParamFamily")) useLast <- TRUE
#            if(useLast && !is(modifyIC(IC), "NULL") ){
#                newParam <- param(L2Fam)
#                main(newParam)[] <- res
#                newL2Fam <- modifyModel(L2Fam, newParam)
#                IC <- modifyIC(IC)(newL2Fam, IC)
#                Infos <- rbind(Infos, c("oneStepEstimator",
#                                        "computation of IC, asVar and asBias via useLast = TRUE"))
#            }else{
#                if(useLast && is(modifyIC(IC), "NULL")){
#                    warning("'useLast = TRUE' only possible if slot 'modifyIC' of 'IC'
#                             is filled with some function!")
#                    Infos <- rbind(Infos, c("oneStepEstimator",
#                                            "slot 'modifyIC' of 'IC' was not filled!"))
#                }
#                Infos <- rbind(Infos, c("oneStepEstimator",
#                                        "computation of IC, asVar and asBias via useLast = FALSE"))
#            }
#            if("asCov" %in% names(Risks(IC)))
#                if(is.matrix(Risks(IC)$asCov) || length(Risks(IC)$asCov) == 1)
#                    asVar <- Risks(IC)$asCov
#                else
#                    asVar <- Risks(IC)$asCov$value
#            else
#                asVar <- getRiskIC(IC, risk = asCov())$asCov$value
#
#            if("asBias" %in% names(Risks(IC))){
#                if(length(Risks(IC)$asBias) == 1)
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias
#                else
#                    asBias <- neighborRadius(IC)*Risks(IC)$asBias$value
#            }else{
#                if(is(IC, "HampIC")){
#                    r <- neighborRadius(IC)
#                    asBias <- r*getRiskIC(IC, risk = asBias(), neighbor = neighbor(IC))$asBias$value
#                }else{
#                    asBias <- NULL
#                }
#            }
#        }else{
#            Infos <- matrix(c("oneStepEstimator", "1-step estimate"), ncol = 2)
#            colnames(Infos) <- c("method", "message")
#            asVar <- NULL
#            asBias <- NULL
#        }
#
#        new("kStepEstimate", name = "1-step estimate", estimate = res,
#            estimate.call = es.call, samplesize = ncol(x), asvar = asVar,
#            asbias = asBias, pIC = IC, steps = 1L, Infos = Infos)
#    })



















