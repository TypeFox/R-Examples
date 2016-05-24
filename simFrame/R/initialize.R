# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------


# model based data control
setMethod("initialize", "DataControl", 
    function(.Object, ...) {
        args <- list(...)
        # use 'rnorm' as default
        if(is.null(args$distribution)) setDistribution(.Object, rnorm)
        callNextMethod()  # call method for superclass (or default)
    })

# sample control
setMethod("initialize", "SampleControl", 
    function(.Object, ...) {
        args <- list(...)
        # use simple random sampling as default
        if(is.null(args$fun)) setFun(.Object, srs)
        callNextMethod()  # call method for superclass (or default)
    })

# two-stage sample control
setMethod("initialize", "TwoStageControl", 
    function(.Object, ...) {
        args <- list(...)
        # use simple random sampling as default
        if(is.null(args$fun)) setFun(.Object, list(srs, srs))
        callNextMethod()  # call method for superclass (or default)
    })

# contamination distributed completely at random (DCAR)
setMethod("initialize", "DCARContControl", 
    function(.Object, ...) {
        args <- list(...)
        # use standard normal distribution as default for contamination data
        if(is.null(args$distribution)) setDistribution(.Object, rnorm)
        callNextMethod()  # call method for superclass (or default)
    })

## insertion of missing values
#setMethod("initialize", "NAControl", 
#    function(.Object, ...) {
#        args <- list(...)
#        # make sure logical indicator whether NAs should be insereted 
#        # into contaminated observations is either TRUE or FALSE
#        intoContamination <- isTRUE(getIntoContamination(.Object))
#        setIntoContamination(.Object, intoContamination)
#        callNextMethod()  # call method for superclass (or default)
#    })
