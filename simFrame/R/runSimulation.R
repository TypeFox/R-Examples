# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## for convenience: construct "SimControl" object and re-call function
setMethod("runSimulation", 
    signature(x = "ANY", setup = "ANY", nrep = "ANY", control = "missing"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        control <- SimControl(contControl=contControl, NAControl=NAControl, 
            design=design, fun=fun, dots=list(...), SAE=SAE)
        if(missing(setup)) {
            if(missing(nrep)) runSimulation(x, control=control)
            else runSimulation(x, nrep=nrep, control=control)
        } else {
            if(missing(nrep)) runSimulation(x, setup, control=control)
            else runSimulation(x, setup, nrep, control)
        }
    })


## design-based simulation
setMethod("runSimulation", 
    signature(x = "data.frame", setup = "VirtualSampleControl", 
        nrep = "missing", control = "SimControl"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        setup <- setup(x, setup)
        runSimulation(x, setup, control=control)
    })

setMethod("runSimulation", 
    signature(x = "data.frame", setup = "SampleSetup", 
        nrep = "missing", control = "SimControl"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        # initializations
        nsam <- length(setup)
        design <- getDesign(control)
        if(nsam == 0) {  # nothing to do
            return(SimResults(design=design, sampleControl=getControl(setup), 
                    control=control))
        }
        contControl <- getContControl(control)
        epsilon <- if(is.null(contControl)) numeric() else getEpsilon(contControl)
        NAControl <- getNAControl(control)
        NArate <- if(is.null(NAControl)) numeric() else getNArate(NAControl)
        # run the simulations
        s <- 1:nsam
        tmp <- lapply(s, designSimulation, x, setup, control)
        # construct results
        getSimResults(tmp, s, epsilon=epsilon, NArate=NArate, design=design, 
            sampleControl=getControl(setup), control=control)
    })


## model-based simulation
setMethod("runSimulation", 
    signature(x = "VirtualDataControl", setup = "missing", 
        nrep = "numeric", control = "SimControl"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        # initializations
        if(length(nrep) == 0) stop("'nrep' must be a non-negative integer")
        else if(length(nrep) > 1) nrep <- nrep[1]
        design <- getDesign(control)
        if(nrep == 0) {  # nothing to do
            return(SimResults(design=design, dataControl=x, 
                    nrep=nrep, control=control))
        }
        contControl <- getContControl(control)
        epsilon <- if(is.null(contControl)) numeric() else getEpsilon(contControl)
        NAControl <- getNAControl(control)
        NArate <- if(is.null(NAControl)) numeric() else getNArate(NAControl)
        # run the simulations
        r <- 1:nrep
        tmp <- lapply(r, modelSimulation, x, control)
        # construct results
        getSimResults(tmp, reps=r, epsilon=epsilon, NArate=NArate, 
            design=design, dataControl=x, control=control)
    })


## mixed simulation designs
setMethod("runSimulation", 
    signature(x = "VirtualDataControl", setup = "VirtualSampleControl", 
        nrep = "numeric", control = "SimControl"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        # initializations
        if(length(nrep) == 0) stop("'nrep' must be a non-negative integer")
        else if(length(nrep) > 1) nrep <- nrep[1]
        nsam <- length(setup)
        design <- getDesign(control)
        if(nrep == 0 || nsam == 0) {  # nothing to do
            return(SimResults(design=design, dataControl=x, 
                    sampleControl=setup, nrep=nrep, control=control))
        }
        contControl <- getContControl(control)
        epsilon <- if(is.null(contControl)) numeric() else getEpsilon(contControl)
        NAControl <- getNAControl(control)
        NArate <- if(is.null(NAControl)) numeric() else getNArate(NAControl)
        # run the simulations (generate data repeatedly and draw samples)
        r <- 1:nrep
        tmp <- lapply(r, mixedSimulation, x, setup, control)
        # construct results
        getSimResults(tmp, 1:nsam, r, epsilon=epsilon, NArate=NArate, 
            design=design, dataControl=x, sampleControl=setup, control=control)
    })


## simulation with repetitions based on (possibly) real data
setMethod("runSimulation",
    signature(x = "data.frame", setup = "missing", 
        nrep = "numeric", control = "SimControl"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        # initializations
        if(length(nrep) == 0) stop("'nrep' must be a non-negative integer")
        else if(length(nrep) > 1) nrep <- nrep[1]
        design <- getDesign(control)
        if(nrep == 0) {  # nothing to do
            return(SimResults(design=design, nrep=nrep, control=control))
        }
        contControl <- getContControl(control)
        epsilon <- if(is.null(contControl)) numeric() else getEpsilon(contControl)
        NAControl <- getNAControl(control)
        NArate <- if(is.null(NAControl)) numeric() else getNArate(NAControl)
        SAE <- isTRUE(getSAE(control))
        # get results (adjustments are needed for parallel computing)
        if(length(design)) {
            spl <- getStrataSplit(x, design, USE.NAMES=FALSE)
            leg <- getStrataLegend(x, design)
            if(SAE) {
                tmp <- replicate(nrep, 
                    manageSimulationSAE(x, spl, control, leg),
                    simplify=FALSE)
            } else {
                xSpl <- lapply(spl, function(s, x) x[s, , drop=FALSE], x)
                tmp <- replicate(nrep, 
                    manageSimulationStrata(xSpl, spl, control, leg), 
                    simplify=FALSE)
            }
        } else {
            tmp <- replicate(nrep, manageSimulation(x, control), 
                simplify=FALSE)
        }
        # construct results
        getSimResults(tmp, reps=1:nrep, epsilon=epsilon, 
            NArate=NArate, design=design, control=control)
    })


## no samples, no repetitions (probably useless, but for completeness)
setMethod("runSimulation",
    signature(x = "data.frame", setup = "missing", 
        nrep = "missing", control = "SimControl"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        runSimulation(x, nrep=1, control=control)
    })

setMethod("runSimulation",
    signature(x = "VirtualDataControl", setup = "missing",
        nrep = "missing", control = "SimControl"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., SAE = FALSE) {
        runSimulation(x, nrep=1, control=control)
    })

setMethod("runSimulation",
    signature(x = "VirtualDataControl", setup = "VirtualSampleControl",
        nrep = "missing", control = "SimControl"),
    function(x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., SAE = FALSE) {
        runSimulation(x, setup, nrep=1, control=control)
    })


## these functions need to be available in the namespace for parallel computing

designSimulation <- function(i, x, setup, control) {
#    print(paste(Sys.time(), i, sep=": "))
    sam <- drawS3(x, setup, i)
    if(nrow(sam) == 0) return(getEmptyResults(control))
    design <- getDesign(control)
    SAE <- isTRUE(getSAE(control))
    if(length(design)) {
        spl <- getStrataSplit(sam, design, USE.NAMES=FALSE)
        leg <- getStrataLegend(sam, design)
        if(SAE) manageSimulationSAE(sam, spl, control, leg)
        else {
            samSpl <- lapply(spl, function(s, x) x[s, , drop=FALSE], sam)
            manageSimulationStrata(samSpl, spl, control, leg)
        }
    } else manageSimulation(sam, control)
}

# 'i' as first argument is necessary for parallel computing  with 'parLapply'
modelSimulation <- function(i, x, control) {
    md <- try(generate(x))
    if(class(md) == "try-error") return(getEmptyResults(control))
    design <- getDesign(control)
    SAE <- isTRUE(getSAE(control))
    if(length(design)) {
        spl <- getStrataSplit(md, design, USE.NAMES=FALSE)
        leg <- getStrataLegend(md, design)
        if(SAE) manageSimulationSAE(md, spl, control, leg)
        else {
            mdSpl <- lapply(spl, function(s, x) x[s, , drop=FALSE], md)
            manageSimulationStrata(mdSpl, spl, control, leg)
        }
    } else manageSimulation(md, control)
}

mixedSimulation <- function(i, x, setup, control) {
    nsam <- length(setup)
    # generate data
    md <- try(generate(x))
    if(class(md) == "try-error") {
        # return empty results for each sample
        return(replicate(nsam, getEmptyResults(control), simplify=FALSE))
    }
    # set up samples (empty vector in case of error in one iteration)
    setup <- setup(md, setup)
    # run design-based simulations on the generated data with the set up samples
    s <- 1:nsam
    lapply(s, designSimulation, md, setup, control)
}

manageSimulation <- function(x, control) {
    # initializations
    contControl <- getContControl(control)
    neps <- length(contControl)
    NAControl <- getNAControl(control)
    nNA <- length(NAControl)
    fun <- getFun(control)
    useOrig <- "orig" %in% argNames(fun)
    dots <- getDots(control)
    # get results
    if(neps) {
        if(nNA) {
            # contamination, missings
            tmp <- lapply(1:neps, 
                function(e) {
                    cx <- try(contaminate(x, contControl, e))
                    if(class(cx) == "try-error") return(vector("list", nNA))
                    lapply(1:nNA, 
                        function(n) try({
                                    nx <- setNA(cx, NAControl, n)
                                    ca <- as.call(c(fun, dots))
                                    ca$x <- nx
                                    if(useOrig) ca$orig <- x
                                    getSimResult(eval(ca))
                                }))
                })
            do.call("c", tmp)
        } else {
            # contamination, no missings
            lapply(1:neps, 
                function(e) try({
                            cx <- contaminate(x, contControl, e)
                            ca <- as.call(c(fun, dots))
                            ca$x <- cx
                            if(useOrig) ca$orig <- x
                            getSimResult(eval(ca))
                        }))
        }
    } else {
        if(nNA) {
            # no contamination, missings
            lapply(1:nNA, 
                function(n) try({
                            nx <- setNA(x, NAControl, n)
                            ca <- as.call(c(fun, dots))
                            ca$x <- nx
                            if(useOrig) ca$orig <- x
                            getSimResult(eval(ca))
                        }))
        } else {
            # no contamination, no missings
            try({
                    ca <- as.call(c(fun, dots))
                    ca$x <- x
                    if(useOrig) ca$orig <- x
                    getSimResult(eval(ca))
                })
        }
    }
}

manageSimulationStrata <- function(xs, indices, control, legend) {
    # initializations
    contControl <- getContControl(control)
    neps <- length(contControl)
    NAControl <- getNAControl(control)
    nNA <- length(NAControl)
    fun <- getFun(control)
    nam <- argNames(fun)
    useOrig <- "orig" %in% nam
    useDomain <- "domain" %in% nam
    dots <- getDots(control)
    # get results
    if(neps) {
        if(nNA) {
            # contamination, missings
            tmp <- lapply(1:neps, 
                function(e) {
                    cxs <- try(lapply(xs, contaminate, contControl, e))
                    if(class(cxs) == "try-error") return(vector("list", nNA))
                    lapply(1:nNA, 
                        function(n) {
                            try({
                                    tmp <- mapply(
                                        function(cx, x, i) {
                                            nx <- setNA(cx, NAControl, n)
                                            ca <- as.call(c(fun, dots))
                                            ca$x <- nx
                                            if(useOrig) ca$orig <- x
                                            if(useDomain) ca$domain <- i
                                            getSimResult(eval(ca))
                                        }, cxs, xs, indices, 
                                        SIMPLIFY=FALSE, USE.NAMES=FALSE)
                                    getSimResultStrata(tmp, legend)
                                })
                        })
                })
            do.call("c", tmp)
        } else {
            # contamination, no missings
            lapply(1:neps, 
                function(e) {
                    try({
                            tmp <- mapply(
                                function(x, i) {
                                    cx <- contaminate(x, contControl, e)
                                    ca <- as.call(c(fun, dots))
                                    ca$x <- cx
                                    if(useOrig) ca$orig <- x
                                    if(useDomain) ca$domain <- i
                                    getSimResult(eval(ca))
                                }, xs, indices, 
                                SIMPLIFY=FALSE, USE.NAMES=FALSE)
                            getSimResultStrata(tmp, legend)
                        })
                })
        }
    } else {
        if(nNA) {
            # no contamination, missings
            lapply(1:nNA, 
                function(n) {
                    try({
                            tmp <- mapply(
                                function(x, i) {
                                    nx <- setNA(x, NAControl, n)
                                    ca <- as.call(c(fun, dots))
                                    ca$x <- nx
                                    if(useOrig) ca$orig <- x
                                    if(useDomain) ca$domain <- i
                                    getSimResult(eval(ca))
                                }, xs, indices, 
                                SIMPLIFY=FALSE, USE.NAMES=FALSE)
                            getSimResultStrata(tmp, legend)
                        })
                })
        } else {
            # no contamination, no missings
            try({
                    tmp <- mapply(
                        function(x, i) {
                            ca <- as.call(c(fun, dots))
                            ca$x <- x
                            if(useOrig) ca$orig <- x
                            if(useDomain) ca$domain <- i
                            getSimResult(eval(ca))
                        }, xs, indices, SIMPLIFY=FALSE, USE.NAMES=FALSE)
                    getSimResultStrata(tmp, legend)
                })
        }
    }
}

manageSimulationSAE <- function(x, indices, control, legend) {
    # initializations
    contControl <- getContControl(control)
    neps <- length(contControl)
    NAControl <- getNAControl(control)
    nNA <- length(NAControl)
    fun <- getFun(control)
    useOrig <- "orig" %in% argNames(fun)
    dots <- getDots(control)
    # get results
    if(neps) {
        if(nNA) {
            # contamination, missings
            tmp <- lapply(1:neps, 
                function(e) {
                    cx <- try(contaminate(x, contControl, e))
                    if(class(cx) == "try-error") return(vector("list", nNA))
                    lapply(1:nNA, 
                        function(n) {
                            try({
                                    nx <- setNA(cx, NAControl, n)
                                    ca <- as.call(c(fun, dots))
                                    ca$x <- nx
                                    if(useOrig) ca$orig <- x
                                    tmp <- lapply(indices, 
                                        function(i) {
                                            ca$domain <- i
                                            getSimResult(eval(ca))
                                        })
                                    getSimResultStrata(tmp, legend)
                                })
                        })
                })
            do.call("c", tmp)
        } else {
            # contamination, no missings
            lapply(1:neps, 
                function(e) {
                    try({
                            cx <- contaminate(x, contControl, e)
                            ca <- as.call(c(fun, dots))
                            ca$x <- cx
                            if(useOrig) ca$orig <- x
                            tmp <- lapply(indices, 
                                function(i) {
                                    ca$domain <- i
                                    getSimResult(eval(ca))
                                })
                            getSimResultStrata(tmp, legend)
                        })
                })
        }
    } else {
        if(nNA) {
            # no contamination, missings
            lapply(1:nNA, 
                function(n) {
                    try({
                            nx <- setNA(x, NAControl, n)
                            ca <- as.call(c(fun, dots))
                            ca$x <- nx
                            if(useOrig) ca$orig <- x
                            tmp <- lapply(indices, 
                                function(i) {
                                    ca$domain <- i
                                    getSimResult(eval(ca))
                                })
                            getSimResultStrata(tmp, legend)
                        })
                })
        } else {
            # no contamination, no missings
            try({
                    ca <- as.call(c(fun, dots))
                    ca$x <- x
                    if(useOrig) ca$orig <- x
                    tmp <- lapply(indices, 
                        function(i) {
                            ca$domain <- i
                            getSimResult(eval(ca))
                        })
                    getSimResultStrata(tmp, legend)
                })
        }
    }
}


## convenience wrapper
runSim <- function(...) {
    res <- runSimulation(...)
    call <- match.call()
    setCall(res, call)
    res
}


## utilities

# check for errors or empty vectors in a list (for simulation results)
checkOK <- function(x) {
    # x ... list
    if(length(x)) {
        sapply(x, function(x) class(x) != "try-error" && length(x$values))
    } else logical()
}

## get empty results
getEmptyResults <- function(control) {
    neps <- length(getContControl(control))
    if(neps == 0) neps <- 1
    nNA <- length(getNAControl(control))
    if(nNA == 0) nNA <- 1
    replicate(neps*nNA, list(values=numeric()))
}

# get result of one simulation run in the correct format
getSimResult <- function(x) {
    if(is(x, "numeric")) x <- list(values=x)
    else if(is(x, "list")) {
        if(!setequal(names(x), c("values", "add"))) {
            stop("the components of the list returned by ", 
                "'fun' must have names 'values' and 'add'")
        }
#    } else if(is(x, "SimResult")) {
#        x <- list(values=getValues(x), add=getAdd(x))
    } else {
#        stop("'fun' must return a numeric vector ", 
#            "or an object of class \"SimResult\"")
        stop("'fun' must return a numeric vector or a ", 
            "list with components 'values' and 'add'")
    }
    x
}

# get result for one stratified simulation run
getSimResultStrata <- function(x, legend) {
    values <- lapply(x, function(x) x$values)
    values <- do.call("rbind", values)
    values <- cbind(legend, values)
    add <- lapply(x, function(x) x$add)
    list(values=values, add=add)
}

# contruct object to be returned
getSimResults <- function(x, samples = numeric(), reps = numeric(), 
    epsilon = numeric(), NArate = numeric(), design = character(), 
    dataControl = NULL, sampleControl = NULL, control) {
    nsam <- length(samples)
    .nrep <- if(missing(reps)) numeric() else length(reps)
    nrep <- length(reps)
    neps <- length(epsilon)
    origNArate <- NArate
    NArate <- convertNArate(NArate)
    nNA <- length(NArate)
    if(is.na(nNA)) nNA <- 0
    # combine results from all runs into one list
    if(nsam && nrep) x <- do.call("c", x)
    if(neps || nNA) x <- do.call("c", x)
    # check for errors or empty results
    nruns <- length(x)
    ok <- checkOK(x)
    x <- x[ok]
    if(length(x) == 0) stop("error or empty result in every simulation run")
    # get additional information (at least one of 'nrep' or 'nsam' is positive)
    ca <- call("expand.grid")  # initialize call
    if(nNA) ca$NArate <- NArate
    if(neps) ca$Epsilon <- epsilon
    if(nsam) ca$Sample <- samples
    if(nrep) ca$Rep <- reps
    info <- eval(ca)  # create data.frame with additional information
    info <- info[, rev(1:ncol(info)), drop=FALSE]  # reverse column order
    info <- cbind(Run=1:nruns, info)[ok, , drop=FALSE]  # add runs
    # elements of 'tmp' are currently lists (faster than S4 objects)
    values <- lapply(x, function(x) x$values)
    if(length(design)) {
        # additional information needs to be adjusted
        reps <- sapply(values, getRepetitions)
#        info <- info[rep(1:nrow(info), each=reps), , drop=FALSE]
        info <- info[rep.int(1:nrow(info), reps), , drop=FALSE]
    }
    values <- do.call("rbind", values)
    # put it all together
    values <- cbind(info, values)
    ninfo <- ncol(info) + length(design)
    nam <- names(values)[-(1:ninfo)]
    rownames(values) <- NULL
    add <- lapply(x, function(x) x$add)
    # return results
    SimResults(values=values, add=add, design=design, colnames=nam, 
        epsilon=epsilon, NArate=origNArate, dataControl=dataControl, 
        sampleControl=sampleControl, nrep=.nrep, control=control)
}
