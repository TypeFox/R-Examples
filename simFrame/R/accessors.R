# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

# mutator functions, whose generic functions contain the '...' argument:
# the expression needs to be defined in the environment of the generic 
# function, but it needs to be evaluated in the environment one level 
# further up (i.e., second environment above the current one)


## class "DataControl"

setMethod("getSize", "DataControl", function(x) slot(x, "size"))
#setMethod("setSize", "DataControl", 
#    function(x, size) eval.parent(substitute(slot(x, "size") <- size)))
setMethod("setSize", "DataControl", 
    function(x, size) {
        eval.parent(substitute(slot(x, "size") <- size, env=parent.frame()), n=2)
    })

setMethod("getDistribution", "DataControl", function(x) slot(x, "distribution"))
setMethod("setDistribution", "DataControl", 
    function(x, distribution) {
        eval.parent(substitute(slot(x, "distribution") <- distribution))
    })

setMethod("getDots", "DataControl", function(x) slot(x, "dots"))
#setMethod("setDots", "DataControl", 
#    function(x, dots) eval.parent(substitute(slot(x, "dots") <- dots)))
setMethod("setDots", "DataControl", 
    function(x, dots) {
        eval.parent(substitute(slot(x, "dots") <- dots, env=parent.frame()), n=2)
    })

setMethod("getColnames", "DataControl", function(x) slot(x, "colnames"))
setMethod("setColnames", "DataControl", 
    function(x, colnames) {
        eval.parent(substitute(slot(x, "colnames") <- colnames))
    })


## class "SampleControl"

setMethod("getK", "VirtualSampleControl", function(x) slot(x, "k"))
setMethod("setK", "VirtualSampleControl", 
    function(x, k) eval.parent(substitute(slot(x, "k") <- k)))

setMethod("getDesign", "SampleControl", function(x) slot(x, "design"))
setMethod("setDesign", "SampleControl", 
    function(x, design) eval.parent(substitute(slot(x, "design") <- design)))

setMethod("getGrouping", "SampleControl", function(x) slot(x, "grouping"))
setMethod("setGrouping", "SampleControl", 
    function(x, grouping) {
        eval.parent(substitute(slot(x, "grouping") <- grouping))
    })

setMethod("getCollect", "SampleControl", function(x) slot(x, "collect"))
setMethod("setCollect", "SampleControl", 
    function(x, collect) eval.parent(substitute(slot(x, "collect") <- collect)))

setMethod("getFun", "SampleControl", function(x) slot(x, "fun"))
#setMethod("setFun", "SampleControl", 
#    function(x, fun) eval.parent(substitute(slot(x, "fun") <- fun)))
setMethod("setFun", "SampleControl", 
    function(x, fun) {
        eval.parent(substitute(slot(x, "fun") <- fun, env=parent.frame()), n=2)
    })

setMethod("getSize", "SampleControl", function(x) slot(x, "size"))
#setMethod("setSize", "SampleControl", 
#    function(x, size) eval.parent(substitute(slot(x, "size") <- size)))
setMethod("setSize", "SampleControl", 
    function(x, size) {
        eval.parent(substitute(slot(x, "size") <- size, env=parent.frame()), n=2)
    })

setMethod("getProb", "SampleControl", function(x) slot(x, "prob"))
#setMethod("setProb", "SampleControl", 
#    function(x, prob) eval.parent(substitute(slot(x, "prob") <- prob)))
setMethod("setProb", "SampleControl", 
    function(x, prob) {
        eval.parent(substitute(slot(x, "prob") <- prob, env=parent.frame()), n=2)
    })

setMethod("getDots", "SampleControl", function(x) slot(x, "dots"))
#setMethod("setDots", "SampleControl", 
#    function(x, dots) eval.parent(substitute(slot(x, "dots") <- dots)))
setMethod("setDots", "SampleControl", 
    function(x, dots) {
        eval.parent(substitute(slot(x, "dots") <- dots, env=parent.frame()), n=2)
    })


## class "TwoStageControl"

setMethod("getDesign", "TwoStageControl", function(x) slot(x, "design"))
setMethod("setDesign", "TwoStageControl", 
    function(x, design) eval.parent(substitute(slot(x, "design") <- design)))

setMethod("getGrouping", "TwoStageControl", function(x) slot(x, "grouping"))
setMethod("setGrouping", "TwoStageControl", 
    function(x, grouping) {
        eval.parent(substitute(slot(x, "grouping") <- grouping))
    })


# utility function to check the 'stage' argument of the following methods
checkStage <- function(stage) {
    if(!isTRUE(stage == 1) && !isTRUE(stage == 2)) {
        stop("'stage' must be either 1 or 2")
    }
}

# in the following mutators: 'stage' is not available in the environment of 
# the generic function and needs to be extracted from the additional arguments

setMethod("getFun", "TwoStageControl", 
    function(x, stage = NULL) {
        fun <- slot(x, "fun")
        if(is.null(stage)) fun
        else {
            checkStage(stage)
            fun[[stage]]
        }
    })
setMethod("setFun", "TwoStageControl", 
    function(x, fun, stage = NULL) {
        pf <- parent.frame()  # environment of generic function
        if(is.null(stage)) expr <- substitute(slot(x, "fun") <- fun, pf)
        else {
            checkStage(stage)
            expr <- substitute(slot(x, "fun")[[list(...)$stage]] <- fun, pf)
        }
        eval.parent(expr, n=2)  # evaluate expression
    })

setMethod("getSize", "TwoStageControl", 
    function(x, stage = NULL) {
        size <- slot(x, "size")
        if(is.null(stage)) size
        else {
            checkStage(stage)
            size[[stage]]
        }
    })
setMethod("setSize", "TwoStageControl", 
    function(x, size, stage = NULL) {
        pf <- parent.frame()  # environment of generic function
        if(is.null(stage)) expr <- substitute(slot(x, "size") <- size, pf)
        else {
            checkStage(stage)
            expr <- substitute(slot(x, "size")[[list(...)$stage]] <- size, pf)
        }
        eval.parent(expr, n=2)  # evaluate expression
    })

setMethod("getProb", "TwoStageControl", 
    function(x, stage = NULL) {
        prob <- slot(x, "prob")
        if(is.null(stage)) prob
        else {
            checkStage(stage)
            prob[[stage]]
        }
    })
setMethod("setProb", "TwoStageControl", 
    function(x, prob, stage = NULL) {
        pf <- parent.frame()  # environment of generic function
        if(is.null(stage)) expr <- substitute(slot(x, "prob") <- prob, pf)
        else {
            checkStage(stage)
            expr <- substitute(slot(x, "prob")[[list(...)$stage]] <- prob, pf)
        }
        eval.parent(expr, n=2)  # evaluate expression
    })

setMethod("getDots", "TwoStageControl", 
    function(x, stage = NULL) {
        dots <- slot(x, "dots")
        if(is.null(stage)) dots
        else {
            checkStage(stage)
            dots[[stage]]
        }
    })
setMethod("setDots", "TwoStageControl", 
    function(x, dots, stage = NULL) {
        pf <- parent.frame()  # environment of generic function
        if(is.null(stage)) expr <- substitute(slot(x, "dots") <- dots, pf)
        else {
            checkStage(stage)
            expr <- substitute(slot(x, "dots")[[list(...)$stage]] <- dots, pf)
        }
        eval.parent(expr, n=2)  # evaluate expression
    })


## class "SampleSetup"

# public accessors (getters)
setMethod("getIndices", "SampleSetup", function(x) slot(x, "indices"))
setMethod("getProb", "SampleSetup", function(x) slot(x, "prob"))
#setMethod("getDesign", "SampleSetup", function(x) slot(x, "design"))
#setMethod("getGrouping", "SampleSetup", function(x) slot(x, "grouping"))
#setMethod("getCollect", "SampleSetup", function(x) slot(x, "collect"))
#setMethod("getFun", "SampleSetup", function(x) slot(x, "fun"))
setMethod("getControl", "SampleSetup", function(x) slot(x, "control"))
setMethod("getSeed", "SampleSetup", function(x) slot(x, "seed"))
setMethod("getCall", "SampleSetup", function(x) slot(x, "call"))

# private mutators (setters)
setMethod("setIndices", "SampleSetup", 
    function(x, indices) eval.parent(substitute(slot(x, "indices") <- indices)))
setMethod("setSeed", "SampleSetup", 
    function(x, seed) eval.parent(substitute(slot(x, "seed") <- seed)))
setMethod("setCall", "SampleSetup", 
    function(x, call) eval.parent(substitute(slot(x, "call") <- call)))

# summary
setMethod("getSize", "SummarySampleSetup", function(x) slot(x, "size"))


## class "ContControl"

setMethod("getTarget", "VirtualContControl", function(x) slot(x, "target"))
setMethod("setTarget", "VirtualContControl", 
    function(x, target) eval.parent(substitute(slot(x, "target") <- target)))

setMethod("getEpsilon", "VirtualContControl", function(x) slot(x, "epsilon"))
setMethod("setEpsilon", "VirtualContControl", 
    function(x, epsilon) eval.parent(substitute(slot(x, "epsilon") <- epsilon)))

setMethod("getGrouping", "ContControl", function(x) slot(x, "grouping"))
setMethod("setGrouping", "ContControl", 
    function(x, grouping) {
        eval.parent(substitute(slot(x, "grouping") <- grouping))
    })

setMethod("getAux", "ContControl", function(x) slot(x, "aux"))
setMethod("setAux", "ContControl", 
    function(x, aux) eval.parent(substitute(slot(x, "aux") <- aux)))

setMethod("getDistribution", "DCARContControl", 
    function(x) slot(x, "distribution"))
setMethod("setDistribution", "DCARContControl", 
    function(x, distribution) {
        eval.parent(substitute(slot(x, "distribution") <- distribution))
    })

setMethod("getDots", "DCARContControl", function(x) slot(x, "dots"))
#setMethod("setDots", "DCARContControl", 
#    function(x, dots) eval.parent(substitute(slot(x, "dots") <- dots)))
setMethod("setDots", "DCARContControl", 
    function(x, dots) {
        eval.parent(substitute(slot(x, "dots") <- dots, env=parent.frame()), n=2)
    })

setMethod("getFun", "DARContControl", function(x) slot(x, "fun"))
#setMethod("setFun", "DARContControl", 
#    function(x, fun) eval.parent(substitute(slot(x, "fun") <- fun)))
setMethod("setFun", "DARContControl", 
    function(x, fun) {
        eval.parent(substitute(slot(x, "fun") <- fun, env=parent.frame()), n=2)
    })

setMethod("getDots", "DARContControl", function(x) slot(x, "dots"))
#setMethod("setDots", "DARContControl", 
#    function(x, dots) eval.parent(substitute(slot(x, "dots") <- dots)))
setMethod("setDots", "DARContControl", 
    function(x, dots) {
        eval.parent(substitute(slot(x, "dots") <- dots, env=parent.frame()), n=2)
    })


## class "NAControl"

setMethod("getTarget", "VirtualNAControl", function(x) slot(x, "target"))
setMethod("setTarget", "VirtualNAControl", 
    function(x, target) eval.parent(substitute(slot(x, "target") <- target)))

setMethod("getNArate", "VirtualNAControl", function(x) slot(x, "NArate"))
setMethod("setNArate", "VirtualNAControl", 
    function(x, NArate) eval.parent(substitute(slot(x, "NArate") <- NArate)))

setMethod("getGrouping", "NAControl", function(x) slot(x, "grouping"))
setMethod("setGrouping", "NAControl", 
    function(x, grouping) {
        eval.parent(substitute(slot(x, "grouping") <- grouping))
    })

setMethod("getAux", "NAControl", function(x) slot(x, "aux"))
setMethod("setAux", "NAControl", 
    function(x, aux) eval.parent(substitute(slot(x, "aux") <- aux)))

setMethod("getIntoContamination", "NAControl", 
    function(x) slot(x, "intoContamination"))
setMethod("setIntoContamination", "NAControl", 
    function(x, intoContamination) {
        eval.parent(substitute(slot(x, "intoContamination") <- intoContamination))
    })


## class "Strata"

# public accessors (getters)
setMethod("getValues", "Strata", function(x) slot(x, "values"))
setMethod("getSplit", "Strata", function(x) slot(x, "split"))
setMethod("getDesign", "Strata", function(x) slot(x, "design"))
setMethod("getNr", "Strata", function(x) slot(x, "nr"))
setMethod("getLegend", "Strata", function(x) slot(x, "legend"))
setMethod("getSize", "Strata", function(x) slot(x, "size"))
setMethod("getCall", "Strata", function(x) slot(x, "call"))

# private mutators (setters)
setMethod("setCall", "Strata", 
    function(x, call) eval.parent(substitute(slot(x, "call") <- call)))


## class "SimControl"

setMethod("getContControl", "SimControl", function(x) slot(x, "contControl"))
setMethod("setContControl", "SimControl", 
    function(x, contControl) {
        eval.parent(substitute(slot(x, "contControl") <- contControl))
    })

setMethod("getNAControl", "SimControl", function(x) slot(x, "NAControl"))
setMethod("setNAControl", "SimControl", 
    function(x, NAControl) {
        eval.parent(substitute(slot(x, "NAControl") <- NAControl))
    })

setMethod("getDesign", "SimControl", function(x) slot(x, "design"))
setMethod("setDesign", "SimControl", 
    function(x, design) eval.parent(substitute(slot(x, "design") <- design)))

setMethod("getFun", "SimControl", function(x) slot(x, "fun"))
#setMethod("setFun", "SimControl", 
#    function(x, fun) eval.parent(substitute(slot(x, "fun") <- fun)))
setMethod("setFun", "SimControl", 
    function(x, fun) {
        eval.parent(substitute(slot(x, "fun") <- fun, env=parent.frame()), n=2)
    })

setMethod("getDots", "SimControl", function(x) slot(x, "dots"))
#setMethod("setDots", "SimControl", 
#    function(x, dots) eval.parent(substitute(slot(x, "dots") <- dots)))
setMethod("setDots", "SimControl", 
    function(x, dots) {
        eval.parent(substitute(slot(x, "dots") <- dots, env=parent.frame()), n=2)
    })

setMethod("getSAE", "SimControl", function(x) slot(x, "SAE"))
setMethod("setSAE", "SimControl", 
    function(x, SAE) eval.parent(substitute(slot(x, "SAE") <- SAE)))


### class "SimResult"
#
## public accessors (getters)
#setMethod("getValues", "SimResult", function(x) slot(x, "values"))
#setMethod("getAdd", "SimResult", function(x) slot(x, "add"))


## class "SimResults"

# public accessors (getters)
setMethod("getValues", "SimResults", function(x) slot(x, "values"))
setMethod("getAdd", "SimResults", function(x) slot(x, "add"))
setMethod("getDesign", "SimResults", function(x) slot(x, "design"))
setMethod("getColnames", "SimResults", function(x) slot(x, "colnames"))
setMethod("getEpsilon", "SimResults", function(x) slot(x, "epsilon"))
setMethod("getNArate", "SimResults", function(x) slot(x, "NArate"))
setMethod("getDataControl", "SimResults", function(x) slot(x, "dataControl"))
setMethod("getSampleControl", "SimResults", function(x) slot(x, "sampleControl"))
setMethod("getNrep", "SimResults", function(x) slot(x, "nrep"))
setMethod("getControl", "SimResults", function(x) slot(x, "control"))
setMethod("getSeed", "SimResults", function(x) slot(x, "seed"))
setMethod("getCall", "SimResults", function(x) slot(x, "call"))

# private mutators (setters)
setMethod("setValues", "SimResults", 
    function(x, values) eval.parent(substitute(slot(x, "values") <- values)))
setMethod("setSeed", "SimResults", 
    function(x, seed) eval.parent(substitute(slot(x, "seed") <- seed)))
setMethod("setCall", "SimResults", 
    function(x, call) eval.parent(substitute(slot(x, "call") <- call)))
