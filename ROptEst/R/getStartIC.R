setMethod("getStartIC",signature(model = "ANY", risk = "ANY"),
           function(model, risk, ...) stop("not yet implemented"))

setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "asGRisk"),
           function(model, risk, ..., withEvalAsVar = TRUE, ..debug=FALSE){
    mc <- match.call(expand.dots=FALSE, call = sys.call(sys.parent(1)))
    dots <- as.list(mc$"...")
    if("fsCor" %in% names(dots)){
        fsCor <- eval(dots[["fsCor"]])
        dots$fsCor <- NULL
    }else fsCor <- 1
    if("eps" %in% names(dots)){
       eps <- dots[["eps"]]
       dots$eps <- NULL
    }else eps <- NULL
    if("neighbor" %in% names(dots)){
       neighbor <- eval(dots[["neighbor"]])
       dots$neighbor <- NULL
    }else neighbor <- ContNeighborhood()

    sm.rmx <- selectMethod("radiusMinimaxIC", signature(
               class(model),class(neighbor),class(risk)))
    dots.rmx <- .fix.in.defaults(dots, sm.rmx)
    dots.rmx$L2Fam <- NULL
    dots.rmx$neighbor <- NULL
    dots.rmx$risk <- NULL

    infMod <- InfRobModel(center = model, neighbor = neighbor)
    sm.optic <- selectMethod("optIC", signature(
                                     class(infMod),class(risk)))
    dots.optic <- .fix.in.defaults(dots, sm.optic)
    dots.optic$model <- NULL
    dots.optic$risk <- NULL
    dots.optic$.withEvalAsVar <- withEvalAsVar

    if(is.null(eps[["e"]])){
        dots.rmx$loRad <- eps$sqn * eps$lower
        dots.rmx$upRad <- eps$sqn * eps$upper
        arg.rmx <- c(list(L2Fam = model, neighbor = neighbor,
                                   risk = risk), dots.rmx)
        if(..debug) print(c(arg.rmx=arg.rmx))
        ICstart <- do.call(radiusMinimaxIC, args=arg.rmx, envir=parent.frame(2))
        if(..debug) print(ICstart)
        if(!isTRUE(all.equal(fsCor, 1, tol = 1e-3))){
            neighbor@radius <- neighborRadius(ICstart)*fsCor
            arg.optic <- c(list( model = infMod, risk = risk), dots.optic)
            if(..debug) print(c(arg.optic=arg.optic))
            ICstart <- do.call(optIC, args = arg.optic, envir=parent.frame(2))
        }
    }else{
        neighbor@radius <- eps$sqn * fsCor * eps$e
        infMod <- InfRobModel(center = model, neighbor = neighbor)
        arg.optic <- c(list(model = infMod, risk = risk),
                           dots.optic)
        if(..debug) print(c(arg.optic=arg.optic))
        ICstart <- do.call(optIC, args = arg.optic, envir=parent.frame(2))
    }
  return(ICstart)
  })

setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "asCov"),
           function(model, risk, ..., ..debug=FALSE){
    return(optIC(model, risk))
  })
setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "trAsCov"),
     getMethod("getStartIC", signature(model = "L2ParamFamily", risk = "asCov"))
           )

setMethod("getStartIC",signature(model = "L2ParamFamily", risk = "asBias"),
           function(model, risk, ..., ..debug=FALSE){
    mc <- match.call(expand.dots=FALSE, call = sys.call(sys.parent(1)))
    dots <- as.list(mc$"...")
    if("neighbor" %in% names(dots)){
       neighbor <- eval(dots[["neighbor"]])
       dots$neighbor <- NULL
    }else neighbor <- ContNeighborhood()

    infMod <- InfRobModel(center = model, neighbor = neighbor)
    return(optIC(infMod, risk))
           })


