# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## for convenience: construct "SimControl" object and re-call function
setMethod("clusterRunSimulation", 
    signature(cl = "ANY", x = "ANY", setup = "ANY", 
        nrep = "ANY", control = "missing"),
    function(cl, x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        control <- SimControl(contControl=contControl, NAControl=NAControl, 
            design=design, fun=fun, dots=list(...), SAE=SAE)
        clusterCall(cl, assign, "control", control, envir=.GlobalEnv)
        if(missing(setup)) {
            if(missing(nrep)) clusterRunSimulation(cl, x, control=control)
            else clusterRunSimulation(cl, x, nrep=nrep, control=control)
        } else {
            if(missing(nrep)) clusterRunSimulation(cl, x, setup, control=control)
            else clusterRunSimulation(cl, x, setup, nrep, control)
        }
    })


## design-based simulation
setMethod("clusterRunSimulation", 
    signature(cl = "ANY", x = "data.frame", setup = "VirtualSampleControl", 
        nrep = "missing", control = "SimControl"),
    function(cl, x, setup, nrep, control, contControl = NULL, 
        NAControl = NULL, design = character(), fun, ..., 
        SAE = FALSE) {
        setup <- clusterSetup(cl, x, setup)
        clusterExport(cl, "setup", envir=sys.frame())
        clusterRunSimulation(cl, x, setup, control=control)
    })

setMethod("clusterRunSimulation", 
    signature(cl = "ANY", x = "data.frame", setup = "SampleSetup", 
        nrep = "missing", control = "SimControl"),
    function(cl, x, setup, nrep, control, contControl = NULL, 
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
        tmp <- parLapply(cl, s, designSimulation, x, setup, control)
        # construct results
        getSimResults(tmp, s, epsilon=epsilon, NArate=NArate, design=design, 
            sampleControl=getControl(setup), control=control)
    })


## model-based simulation
setMethod("clusterRunSimulation", 
    signature(cl = "ANY", x = "VirtualDataControl", 
        setup = "missing", nrep = "numeric", control = "SimControl"),
    function(cl, x, setup, nrep, control, contControl = NULL, 
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
        tmp <- parLapply(cl, r, modelSimulation, x, control)
        # construct results
        getSimResults(tmp, reps=r, epsilon=epsilon, NArate=NArate, 
            design=design, dataControl=x, control=control)
    })


## mixed simulation designs
setMethod("clusterRunSimulation", 
    signature(cl = "ANY", x = "VirtualDataControl", 
        setup = "VirtualSampleControl", nrep = "numeric", 
        control = "SimControl"),
    function(cl, x, setup, nrep, control, contControl = NULL, 
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
        tmp <- parLapply(cl, r, mixedSimulation, x, setup, control)
        # construct results
        getSimResults(tmp, 1:nsam, r, epsilon=epsilon, NArate=NArate, 
            design=design, dataControl=x, sampleControl=setup, control=control)
    })


## simulation with replications based on (possibly) real data
setMethod("clusterRunSimulation",
    signature(cl = "ANY", x = "data.frame", setup = "missing", 
        nrep = "numeric", control = "SimControl"),
    function(cl, x, setup, nrep, control, contControl = NULL, 
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
        # get results
        r <- 1:nrep
        if(length(design)) {
            # necessary objects need to be constructed on workers
            seqList <- clusterSplit(cl, r)
            nrList <- lapply(seqList, length)
            if(SAE) {
                tmp <- clusterApply(cl, nrList, 
                    function(nr, x, control) {
                        design <- getDesign(control)
                        spl <- getStrataSplit(x, design, USE.NAMES=FALSE)
                        leg <- getStrataLegend(x, design)
                        replicate(nr, 
                            manageSimulationSAE(x, spl, control, leg), 
                            simplify=FALSE)
                    }, x, control)
            } else {
                tmp <- clusterApply(cl, nrList, 
                    function(nr, x, control) {
                        design <- getDesign(control)
                        spl <- getStrataSplit(x, design, USE.NAMES=FALSE)
                        xSpl <- lapply(spl, function(s, x) x[s, , drop=FALSE], x)
                        leg <- getStrataLegend(x, design)
                        replicate(nr, 
                            manageSimulationStrata(xSpl, spl, control, leg), 
                            simplify=FALSE)
                    }, x, control)
            }
            tmp <- do.call("c", tmp)
        } else {
            tmp <- parLapply(cl, r, 
                function(i) manageSimulation(x, control))
        }
        # construct results
        getSimResults(tmp, reps=r, epsilon=epsilon, 
            NArate=NArate, design=design, control=control)
    })
