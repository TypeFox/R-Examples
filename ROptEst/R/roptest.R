###############################################################################
## Optimally robust estimation
###############################################################################
roptest.old <- function(x, L2Fam, eps, eps.lower, eps.upper, fsCor = 1, initial.est,
                    neighbor = ContNeighborhood(), risk = asMSE(), steps = 1L, 
                    distance = CvMDist, startPar = NULL, verbose = NULL,
                    OptOrIter = "iterate",
                    useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    na.rm = TRUE, initial.est.ArgList, ...,
                    withLogScale = TRUE){
    if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
    es.call <- match.call()
    if(missing(x))
        stop("'x' is missing with no default")
    if(missing(L2Fam))
        stop("'L2Fam' is missing with no default")
    if(!is.numeric(x)){
        if(is.data.frame(x))
            x <- data.matrix(x)
        else
            x <- as.matrix(x)
        if(!is.matrix(x))
            stop("'x' has to be a numeric vector resp. a matrix or data.frame")
    }
    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)

    if(missing(eps) && missing(eps.lower) && missing(eps.upper)){
        eps.lower <- 0
        eps.upper <- 0.5
    }
    if(missing(eps)){
        if(!missing(eps.lower) && missing(eps.upper))
            eps.upper <- 0.5
        if(missing(eps.lower) && !missing(eps.upper))
            eps.lower <- 0
        if(length(eps.lower) != 1 || length(eps.upper) != 1)
            stop("'eps.lower' and 'eps.upper' have to be of length 1")
        if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper) 
            stop("'eps.lower' < 'eps.upper' is not fulfilled")
        if((eps.lower < 0) || (eps.upper > 0.5))
            stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
    }else{
        if(length(eps) != 1)
            stop("'eps' has to be of length 1")
        if(eps == 0)
            stop("'eps = 0'! => use functions 'mean' and 'sd' for estimation")
        if((eps < 0) || (eps > 0.5))
            stop("'eps' has to be in (0, 0.5]")
    }
    if(fsCor <= 0)
        stop("'fsCor' has to be positive")
    if(length(fsCor) != 1){
        stop("'fsCor' has to be of length 1")
    }
    if(!is.integer(steps))
        steps <- as.integer(steps)
    if(steps < 1){
        stop("'steps' has to be some positive integer value")
    }
    if(length(steps) != 1){
        stop("'steps' has to be of length 1")
    }

    if(missing(initial.est))
        initial.est <- MDEstimator(x = x, ParamFamily = L2Fam, distance = distance,
                                            startPar = startPar, ...)
    nrvalues <-  length(L2Fam@param)
    initial.est <- kStepEstimator.start(initial.est, x = x,
                                        nrvalues = nrvalues, na.rm = na.rm,
                                        L2Fam = L2Fam,
                                        startList = initial.est.ArgList)

    newParam <- param(L2Fam)
    main(newParam)[] <- as.numeric(initial.est)
    L2FamStart <- modifyModel(L2Fam, newParam)
    if(is.matrix(x))
        sqrtn <- sqrt(ncol(x))
    else
        sqrtn <- sqrt(length(x))
    if(missing(eps)){
        r.lower <- sqrtn*eps.lower
        r.upper <- sqrtn*eps.upper
        ICstart <- radiusMinimaxIC(L2Fam = L2FamStart, neighbor = neighbor, risk = risk, 
                                   loRad = r.lower, upRad = r.upper, verbose = verbose,
                                   OptOrIter = OptOrIter, ...)
        if(!isTRUE(all.equal(fsCor, 1, tol = 1e-3))){
            neighbor@radius <- neighborRadius(ICstart)*fsCor
            infMod <- InfRobModel(center = L2FamStart, neighbor = neighbor)
            ICstart <- optIC(model = infMod, risk = risk, verbose = verbose,
                             OptOrIter = OptOrIter, ...)
        }    
    }else{
        neighbor@radius <- sqrtn*eps*fsCor
        infMod <- InfRobModel(center = L2FamStart, neighbor = neighbor)
        ICstart <- optIC(model = infMod, risk = risk, verbose = verbose,
                         OptOrIter = OptOrIter, ...)
    }
    res <- kStepEstimator(x, IC = ICstart, start = initial.est, steps = steps, useLast = useLast,
                          withUpdateInKer = withUpdateInKer, IC.UpdateInKer = IC.UpdateInKer,
                          withICList = withICList, withPICList = withPICList, na.rm = na.rm,
                          withLogScale = withLogScale)
    res@estimate.call <- es.call
    Infos <- matrix(c("roptest", 
                      paste(steps, "-step estimate for ", name(L2Fam), sep = "")),
                    ncol = 2)
    colnames(Infos) <- c("method", "message")

    if(! .isUnitMatrix(trafo(L2Fam)))
       Infos <- rbind(Infos, c("roptest",
                            paste("computation of IC",
                                   ifelse(withUpdateInKer,"with","without") ,
                                   "modification in ker(trafo)")))

    Infos <- rbind(Infos, c("roptest",
                            paste("computation of IC, asvar and asbias via useLast =", useLast)))
    Infos(res) <- Infos
    res@completecases <- completecases
    res@start <- initial.est
    return(res)
}
