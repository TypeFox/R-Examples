# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

setGeneric("clusterRunSimulation",
    function(cl, x, setup, nrep, control, contControl = NULL, 
            NAControl = NULL, design = character(), fun, ..., 
            SAE = FALSE) {
        res <- standardGeneric("clusterRunSimulation")
        call <- match.call()
        setCall(res, call)
        res
    },
    valueClass = "SimResults")

setGeneric("clusterSetup",
    function(cl, x, control, ...) {
        res <- standardGeneric("clusterSetup")
        call <- match.call()
        setCall(res, call)
        res
    },
    valueClass = "SampleSetup")

setGeneric("contaminate",
    function(x, control, ...) standardGeneric("contaminate"),
    valueClass = "data.frame")

setGeneric("draw",
    function(x, setup, ...) standardGeneric("draw"), 
    valueClass = "data.frame")

setGeneric("generate",
    function(control, ...) standardGeneric("generate"), 
    valueClass = "data.frame")

setGeneric("getSampleIndices",
    function(x, control) standardGeneric("getSampleIndices"),
    valueClass = "list")

setGeneric("getSampleProb",
    function(x, control) standardGeneric("getSampleProb"),
    valueClass = "numeric")

setGeneric("getStrataLegend",
    function(x, design) standardGeneric("getStrataLegend"),
    valueClass = "data.frame")

setGeneric("getStrataSplit",
    function(x, design, USE.NAMES = TRUE) standardGeneric("getStrataSplit"),
    valueClass = "list")

setGeneric("getStrataTable",
    function(x, design) standardGeneric("getStrataTable"),
    valueClass = "data.frame")

setGeneric("getStratumSizes",
    function(x, design, USE.NAMES = TRUE) standardGeneric("getStratumSizes"),
    valueClass = "numeric")

setGeneric("getStratumValues",
    function(x, design, split) standardGeneric("getStratumValues"),
    valueClass = "numeric")

setGeneric("runSimulation",
    function(x, setup, nrep, control, contControl = NULL, 
            NAControl = NULL, design = character(), fun, ..., 
            SAE = FALSE) {
        # make sure that .Random.seed exists
        if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
        # call method and store seed before and after
        firstSeed <- .Random.seed
        res <- standardGeneric("runSimulation")
        lastSeed <- .Random.seed
        setSeed(res, list(firstSeed, lastSeed))
        call <- match.call()
        setCall(res, call)
        res
    },
    valueClass = "SimResults")

setGeneric("setNA",
    function(x, control, ...) standardGeneric("setNA"),
    valueClass = "data.frame")

setGeneric("setup",
    function(x, control, ...) {
        # make sure that .Random.seed exists
        if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
        # call method and store seed before and after
        firstSeed <- .Random.seed
        res <- standardGeneric("setup")
        lastSeed <- .Random.seed
        setSeed(res, list(firstSeed, lastSeed))
        call <- match.call()
        setCall(res, call)
        res
    },
    valueClass = "SampleSetup")

setGeneric("simApply",
    function(x, design, fun, ...) standardGeneric("simApply"))

setGeneric("simBwplot",
    function(x, ...) standardGeneric("simBwplot"))

setGeneric("simDensityplot",
    function(x, ...) standardGeneric("simDensityplot"))

setGeneric("simSapply",
    function(x, design, fun, ..., simplify = TRUE) standardGeneric("simSapply"))

setGeneric("simXyplot",
    function(x, ...) standardGeneric("simXyplot"))

setGeneric("stratify", 
    function(x, design) {
        res <- standardGeneric("stratify")
        call <- match.call()
        setCall(res, call)
        res
    }, 
    valueClass = "Strata")


## public accessor and mutator functions (to be exported)
setGeneric("getAdd", function(x) standardGeneric("getAdd"))

setGeneric("getAux", function(x) standardGeneric("getAux"))
setGeneric("setAux", function(x, aux) standardGeneric("setAux"))

setGeneric("getCall", function(x, ...) standardGeneric("getCall"))

setGeneric("getCollect", function(x) standardGeneric("getCollect"))
setGeneric("setCollect", function(x, collect) standardGeneric("setCollect"))

setGeneric("getColnames", function(x) standardGeneric("getColnames"))
setGeneric("setColnames", 
    function(x, colnames) standardGeneric("setColnames"))

setGeneric("getContControl", function(x) standardGeneric("getContControl"))
setGeneric("setContControl", 
    function(x, contControl) standardGeneric("setContControl"))

setGeneric("getControl", function(x) standardGeneric("getControl"))

setGeneric("getDataControl", function(x) standardGeneric("getDataControl"))

setGeneric("getDesign", function(x) standardGeneric("getDesign"))
setGeneric("setDesign", function(x, design) standardGeneric("setDesign"))

setGeneric("getDistribution", function(x) standardGeneric("getDistribution"))
setGeneric("setDistribution", 
    function(x, distribution) standardGeneric("setDistribution"))

#setGeneric("getDots", function(x) standardGeneric("getDots"))
#setGeneric("setDots", function(x, dots) standardGeneric("setDots"))
setGeneric("getDots", function(x, ...) standardGeneric("getDots"))
setGeneric("setDots", function(x, dots, ...) standardGeneric("setDots"))

setGeneric("getEpsilon", function(x) standardGeneric("getEpsilon"))
setGeneric("setEpsilon", function(x, epsilon) standardGeneric("setEpsilon"))

setGeneric("getGrouping", function(x) standardGeneric("getGrouping"))
setGeneric("setGrouping", function(x, grouping) standardGeneric("setGrouping"))

#setGeneric("getFun", function(x) standardGeneric("getFun"))
#setGeneric("setFun", function(x, fun) standardGeneric("setFun"))
setGeneric("getFun", function(x, ...) standardGeneric("getFun"))
setGeneric("setFun", function(x, fun, ...) standardGeneric("setFun"))

setGeneric("getIndices", function(x) standardGeneric("getIndices"))

setGeneric("getIntoContamination", 
    function(x) standardGeneric("getIntoContamination"))
setGeneric("setIntoContamination", 
    function(x, intoContamination) standardGeneric("setIntoContamination"))

setGeneric("getK", function(x) standardGeneric("getK"))
setGeneric("setK", function(x, k) standardGeneric("setK"))

setGeneric("getLegend", function(x) standardGeneric("getLegend"))

setGeneric("getNAControl", function(x) standardGeneric("getNAControl"))
setGeneric("setNAControl", 
    function(x, NAControl) standardGeneric("setNAControl"))

setGeneric("getNArate", function(x) standardGeneric("getNArate"))
setGeneric("setNArate", function(x, NArate) standardGeneric("setNArate"))

setGeneric("getNr", function(x) standardGeneric("getNr"))

setGeneric("getNrep", function(x) standardGeneric("getNrep"))

#setGeneric("getProb", function(x) standardGeneric("getProb"))
#setGeneric("setProb", function(x, prob) standardGeneric("setProb"))
setGeneric("getProb", function(x, ...) standardGeneric("getProb"))
setGeneric("setProb", function(x, prob, ...) standardGeneric("setProb"))

setGeneric("getSAE", function(x) standardGeneric("getSAE"))
setGeneric("setSAE", function(x, SAE) standardGeneric("setSAE"))

setGeneric("getSampleControl", function(x) standardGeneric("getSampleControl"))

setGeneric("getSeed", function(x) standardGeneric("getSeed"))

#setGeneric("getSize", function(x) standardGeneric("getSize"))
#setGeneric("setSize", function(x, size) standardGeneric("setSize"))
setGeneric("getSize", function(x, ...) standardGeneric("getSize"))
setGeneric("setSize", function(x, size, ...) standardGeneric("setSize"))

setGeneric("getSplit", function(x) standardGeneric("getSplit"))

setGeneric("getTarget", function(x) standardGeneric("getTarget"))
setGeneric("setTarget", function(x, target) standardGeneric("setTarget"))

setGeneric("getValues", function(x) standardGeneric("getValues"))


## private accessor and mutator functions (not exported)
setGeneric("setCall", function(x, call) standardGeneric("setCall"))
setGeneric("setIndices", function(x, indices) standardGeneric("setIndices"))
setGeneric("setSeed", function(x, seed) standardGeneric("setSeed"))
setGeneric("setValues", function(x, values) standardGeneric("setValues"))


## existing S3 or S4 generics (just to be safe)
setGeneric("aggregate")
setGeneric("head")
setGeneric("length")
setGeneric("plot")
setGeneric("show")
setGeneric("summary")
setGeneric("tail")
