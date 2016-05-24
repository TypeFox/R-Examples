###############################################################################
## Optimally robust estimation
###############################################################################
.dynScopeEval <- function(expr){
   le <- length(sys.calls())
   i <- 1
   while(i< le){
      a <- try(eval(expr,envir=sys.frame(-i)),silent=TRUE)
      if(!is(a,"try-error")) return(a)
      i <- i + 1
   }
   stop("Could not evaluate expression.")
}

#####################################################################
### probably a dead end: filling in unevaluated arguments is too tricky...
if(FALSE){
roptest <- function(x, L2Fam, eps, eps.lower, eps.upper, fsCor = 1, initial.est,
                    neighbor = ContNeighborhood(), risk = asMSE(), steps = 1L,
                    distance = CvMDist, startPar = NULL, verbose = NULL,
                    OptOrIter = "iterate",
                    useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ...,
                    withLogScale = TRUE,..withCheck=FALSE,
                    withTimings = FALSE, withMDE = NULL,
                    withEvalAsVar = NULL){
    es.call <- es.call.e <- match.call()
    es.call.e <- (as.list(es.call.e))
    es.call.e[["..."]] <- NULL
    for(i in seq(along.with=es.call.e))
        es.call.e[[i]] <- .dynScopeEval(es.call.e[[i]])
    es.call0 <- match.call(expand.dots=FALSE)
    mwt <- !is.null(es.call$withTimings)
    es.call$withTimings <- NULL
    dots <- es.call0[["..."]]
    es.call1 <- .constructArg.list(roptest,es.call.e, onlyFormal=FALSE,
                            debug = ..withCheck)$mc

    res <- .constructArg.list(gennbCtrl,es.call1, onlyFormal=TRUE,
                    debug = ..withCheck)
     nbCtrl0 <-  as.list(res$mc[-1])
     es.call1 <- res$esc
    res <- .constructArg.list(genstartCtrl,es.call1, onlyFormal=TRUE,
                   debug = ..withCheck)
     startCtrl0 <-  as.list(res$mc[-1])
     es.call1 <- res$esc
    res <- .constructArg.list(genkStepCtrl,es.call1, onlyFormal=TRUE,
                   debug = ..withCheck)
     kStepCtrl0 <-  as.list(res$mc[-1])
     es.call1 <- res$esc

    list1 <- as.list(es.call1[-1])
    res <- do.call(gennbCtrl,nbCtrl0)
      if(length(res)) list1$nbCtrl <- res
    res <- do.call(genstartCtrl,startCtrl0)
      if(length(res)) list1$startCtrl <- res
    res <- do.call(genkStepCtrl,kStepCtrl0)
      if(length(res)) list1$kStepCtrl <- res
    list1 <- c(list1,dots)
    list1$..withCheck <- NULL
    list1$withTimings <- withTimings
    if(..withCheck) print(list1)
    if(..withCheck) return(substitute(do.call(robest, lis), list(lis=list1)))
    else{
    res <- do.call(robest, list1)
    if(mwt) es.call$withTimings <- withTimings
    res@estimate.call <- es.call
    return(res)}
}
##############
}
### end of dead end.
#####################################################################

roptest <- function(x, L2Fam, eps, eps.lower, eps.upper, fsCor = 1, initial.est,
                    neighbor = ContNeighborhood(), risk = asMSE(), steps = 1L,
                    distance = CvMDist, startPar = NULL, verbose = NULL,
                    OptOrIter = "iterate",
                    useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ...,
                    withLogScale = TRUE,..withCheck=FALSE,
                    withTimings = FALSE, withMDE = NULL,
                    withEvalAsVar = NULL){
    dots <- match.call(expand.dots=FALSE)[["..."]]
    scalename <- dots[["scalename"]]
    nbCtrl <- list()
    nbCtrl[["neighbor"]] <- if(!missing(neighbor)) neighbor else ContNeighborhood()
    if(!missing(eps)) nbCtrl[["eps"]] <- eps
    if(!missing(eps.lower)) nbCtrl[["eps.lower"]] <- eps.lower
    if(!missing(eps.upper)) nbCtrl[["eps.upper"]] <- eps.upper

    startCtrl <- list()
    if(!missing(initial.est)) startCtrl[["initial.est"]] <- initial.est
    if(!missing(initial.est.ArgList))
       startCtrl[["initial.est.ArgList"]] <- initial.est.ArgList
    startCtrl[["startPar"]] <- if(!missing(startPar)) startPar else NULL
    startCtrl[["distance"]] <- if(!missing(distance)) distance else NULL
    startCtrl[["withMDE"]] <- if(!missing(withMDE)) withMDE else NULL

    kStepCtrl <- list()
    kStepCtrl[["useLast"]] <- if(!missing(useLast)) useLast else getRobAStBaseOption("kStepUseLast")
    kStepCtrl[["withUpdateInKer"]] <- if(!missing(withUpdateInKer)) withUpdateInKer else getRobAStBaseOption("withUpdateInKer")
    kStepCtrl[["IC.UpdateInKer"]] <- if(!missing(IC.UpdateInKer)) IC.UpdateInKer else getRobAStBaseOption("IC.UpdateInKer")
    kStepCtrl[["withICList"]] <- if(!missing(withICList)) withICList else getRobAStBaseOption("withICList")
    kStepCtrl[["withPICList"]] <- if(!missing(withPICList)) withPICList else getRobAStBaseOption("withPICList")
    kStepCtrl[["scalename"]] <- if(!is.null(scalename)) scalename else "scale"
    kStepCtrl[["withLogScale"]] <- if(!missing(withLogScale)) withLogScale else TRUE
    kStepCtrl[["withEvalAsVar"]] <- if(!missing(withEvalAsVar)) withEvalAsVar else NULL

    return(robest(x=x, L2Fam=L2Fam,  fsCor = fsCor,
           risk = risk, steps = steps, verbose = verbose,
           OptOrIter = OptOrIter, nbCtrl = nbCtrl,
           startCtrl = startCtrl, kStepCtrl = kStepCtrl,
           na.rm = na.rm, ..., debug = ..withCheck,
           withTimings = withTimings))
}
#roptest(x=1:10,L2Fam=GammaFamily(),also=3,..withCheck=TRUE)

#roptest(x=1:10,L2Fam=GammaFamily(),also=3)

robest <- function(x, L2Fam,  fsCor = 1,
                     risk = asMSE(), steps = 1L,
                      verbose = NULL,
                    OptOrIter = "iterate",
                    nbCtrl = gennbCtrl(),
                    startCtrl = genstartCtrl(),
                    kStepCtrl = genkStepCtrl(),
                    na.rm = TRUE, ..., debug = FALSE,
                    withTimings = FALSE){

    es.call <- match.call()
    es.call0 <- match.call(expand.dots=FALSE)
    mwt <- !is.null(es.call$withTimings)
    es.call$withTimings <- NULL
    es.call0$withTimings <- NULL
    dots <- es.call0[["..."]]
    es.call0$"..." <- NULL
#    if(!is.null(dots[["escall"]])) es.call <- dots[["escall"]]

    if(missing(nbCtrl)||is.null(nbCtrl))          nbCtrl <- gennbCtrl()
    if(missing(startCtrl)||is.null(startCtrl)) startCtrl <- genstartCtrl()
    if(missing(kStepCtrl)||is.null(kStepCtrl)) kStepCtrl <- genkStepCtrl()
    nbCtrl    <- .fix.in.defaults(nbCtrl, gennbCtrl)
    startCtrl <- .fix.in.defaults(startCtrl, genstartCtrl)
    kStepCtrl <- .fix.in.defaults(kStepCtrl, genkStepCtrl)

    if(missing(L2Fam))
        stop("'L2Fam' is missing with no default")
    if(!is(L2Fam, "L2ParamFamily"))
        stop("'L2Fam' must be of class 'L2ParamFamily'.")

    withEvalAsVar <- kStepCtrl$withEvalAsVar
    if(is.null(withEvalAsVar)) withEvalAsVar <- L2Fam@.withEvalAsVar


    es.list <- as.list(es.call0[-1])
    es.list <- c(es.list,nbCtrl)
    es.list$nbCtrl <- NULL
    es.list0 <-  c(es.list,dots)

    if(missing(verbose)|| is.null(verbose))
           es.list$verbose <- getRobAStBaseOption("all.verbose")

    res.x <- .pretreat(x,na.rm)
    x <- res.x$x
    completecases <- res.x$completecases

    es.list$x <- x
    if(debug){   cat("\nes.list:\n");print(es.list);cat("\n\n")}

    .isOKfsCor(fsCor)

    .isOKsteps(steps)

    if(debug){
      if(is.null(startCtrl$initial.est)){
       print(substitute(MDEstimator(x = x0, ParamFamily = L2Fam0,
                        distance = distance0,
                        startPar = startPar0,dots0),
                        list(x0=x,L2Fam0=L2Fam, distance0=startCtrl$distance,
                             startPar0=startCtrl$startPar0, dots0=dots)))
        startCtrl$initial.est <- "BLUB"
      }
    }else{
      if(is.null(startCtrl$initial.est)){
         startPar0 <- if(is.null(startCtrl$startPar))
                         L2Fam@startPar else startCtrl$startPar
         wMDE <- if(is.null(startCtrl$withMDE))
                         L2Fam@.withMDE else startCtrl$withMDE
         if(is(startPar0, "function")) if(!wMDE){
            startCtrl$initial.est <- function(x,...)startPar0(x)
         }else
            startCtrl$initial.est <- MDEstimator(x = x, ParamFamily = L2Fam,
                                  distance = startCtrl$distance,
                                  startPar = startCtrl$startPar, ...)

      }
    }
    nrvalues <-  length(L2Fam@param)

    if(debug){
      if(!is.null(startCtrl$initial.est.ArgList)){
       initial.est <- print(substitute({kStepEstimator.start(initial.est0, x = x0,
                                        nrvalues = nrvalues0, na.rm = na.rm0,
                                        L2Fam = L2Fam0,
                                        startList = startCtrl0)},
                                        list(x0=x,L2Fam0=L2Fam,
                                        initial.est0=startCtrl$initial.est,
                                        nrvalues0=nrvalues,na.rm0=na.rm,
                                     startCtrl0 = startCtrl$initial.est.ArgList)
                                 ))
      }else{
       print(substitute(kStepEstimator.start(initial.est0, x = x0,
                                        nrvalues = nrvalues0, na.rm = na.rm0,
                                        L2Fam = L2Fam0),list(x0=x,L2Fam0=L2Fam,
                                        initial.est0=startCtrl$initial.est,
                                        nrvalues0=nrvalues,na.rm0=na.rm)
                                 ))
      }
    }else{
    sy.start <- system.time({
      sctrl.init <- eval(startCtrl$initial.est)
      if(is.null(startCtrl$initial.est.ArgList)){
       initial.est <-  kStepEstimator.start(start = sctrl.init, x = x,
                                        nrvalues = nrvalues, na.rm = na.rm,
                                        L2Fam = L2Fam, startList = NULL)
      }else{
       initial.est <-  kStepEstimator.start(start = sctrl.init, x = x,
                                        nrvalues = nrvalues, na.rm = na.rm,
                                        L2Fam = L2Fam,
                                        startList = startCtrl$initial.est.ArgList)
      }
     })
     if(withTimings) print(sy.start)
    }

    newParam <- param(L2Fam)

    if(!debug){
       main(newParam)[] <- as.numeric(initial.est)
       L2FamStart <- modifyModel(L2Fam, newParam)
    }
    if(debug) print(risk)

    if(!is(risk,"interpolRisk"))
       es.list0$eps <- do.call(.check.eps, args=c(nbCtrl,list("x"=x)))
   
    es.list0$risk <- NULL
    es.list0$L2Fam <- NULL
    neighbor <- eval(es.list0$neighbor)
    es.list0$neighbor <- NULL
    es.list0$kStepCtrl <- NULL
    es.list0$startCtrl <- NULL
    es.list0$distance <- NULL
    es.list0$x <- NULL
    es.list0$steps <- NULL
    es.list0$eps.lower <- NULL
    es.list0$eps.upper <- NULL
    es.list0$withEvalAsVar <- NULL
    es.list0$na.rm <- NULL
    es.list0$fsCor <- eval(es.list0$fsCor)

    if(debug) {cat("\n\n\n::::\n\n")
    argList <- c(list(model=L2Fam,risk=risk,neighbor=neighbor),
                                             es.list0)
    print(argList)
    cat("\n\n\n")
    }
    if(!debug){
      sy.getStartIC <-  system.time({
       ICstart <- do.call(getStartIC, args=c(list(model=L2FamStart,risk=risk,
                              neighbor=neighbor, withEvalAsVar = withEvalAsVar),
                              es.list0))
     })
     if (withTimings) print(sy.getStartIC)
     }


      if(debug){
         ICstart <- "BUL"
         argList <- list(x, IC = ICstart, start = initial.est, steps = steps,
                            useLast = kStepCtrl$useLast,
                            withUpdateInKer = kStepCtrl$withUpdateInKer,
                            IC.UpdateInKer = kStepCtrl$IC.UpdateInKer,
                            withICList = kStepCtrl$withICList,
                            withPICList = kStepCtrl$withPICList,
                            na.rm = na.rm,
                            scalename = kStepCtrl$scalename,
                            withLogScale = kStepCtrl$withLogScale,
                            withEvalAsVar = withEvalAsVar)
         print(argList) }
      sy.kStep <- system.time({
         res <- kStepEstimator(x, IC = ICstart, start = initial.est, steps = steps,
                            useLast = kStepCtrl$useLast,
                            withUpdateInKer = kStepCtrl$withUpdateInKer,
                            IC.UpdateInKer = kStepCtrl$IC.UpdateInKer,
                            withICList = kStepCtrl$withICList,
                            withPICList = kStepCtrl$withPICList,
                            na.rm = na.rm,
                            scalename = kStepCtrl$scalename,
                            withLogScale = kStepCtrl$withLogScale,
                            withEvalAsVar = withEvalAsVar)
                            })
       if (withTimings) print(sy.kStep)

    if(!debug){
         if(mwt) es.call$withTimings <- withTimings
         res@estimate.call <- es.call
    }
    Infos <- matrix(c("roptest", 
                      paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                     ncol = 2)
    colnames(Infos) <- c("method", "message")


    if(! .isUnitMatrix(trafo(L2Fam)))
       Infos <- rbind(Infos, c("method"="roptest",
                    "message"=paste("computation of IC",
                               ifelse(kStepCtrl$withUpdateInKer,"with","without") ,
                                   "modification in ker(trafo)")))

    Infos <- rbind(Infos, c("method"="roptest",
     "message"=paste("computation of IC, asvar and asbias via useLast =",
                  kStepCtrl$useLast)))

    if(debug) print(Infos)
    if(debug) return(NULL)
    sy.all <- sy.start+sy.getStartIC+sy.kStep
    if(withTimings) print(sy.all)
    sy <- rbind(start=sy.start,getStartIC=sy.getStartIC,
                kStep=sy.kStep, all=sy.all)
    Infos(res) <- Infos
    res@completecases <- completecases
    res@start <- initial.est
    attr(res, "timings") <- sy
    return(res)
}
