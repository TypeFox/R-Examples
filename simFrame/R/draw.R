# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

setMethod("draw",  
    signature(x = "data.frame", setup = "SampleSetup"), 
    function(x, setup, i = 1) drawS3(x, setup, i))

setMethod("draw", 
    signature(x = "data.frame", setup = "VirtualSampleControl"), 
    function(x, setup) {
        setK(setup, 1)
        draw(x, setup(x, setup), i=1)
    })

setMethod("draw", 
    signature(x = "data.frame", setup = "character"), 
    function(x, setup, ...) {
        if(length(setup) != 1) {
            stop("'setup' must specify exactly one ", 
                "class extending \"VirtualSampleControl\"")
        }
        if(!extends(setup, "VirtualSampleControl")) {
            stop(gettextf("\"%s\" does not extend \"VirtualSampleControl\"", 
                    setup))
        }
#        draw(x, new(setup, ...))
        # temporary solution: constructor for class "TwoStageControl" has 
        # arguments that aren't slots
        if(isTRUE(setup == "TwoStageControl")) setup <- TwoStageControl(...)
        else setup <- new(setup, ...)
        draw(x, setup)
    })

setMethod("draw", 
    signature(x = "data.frame", setup = "missing"), 
    function(x, setup, ...) {
        draw(x, SampleControl(...))
    })


## internal S3 function 
# this is used in 'runSimulation' and 'clusterRunSimulation': there the 
# objects are already checked for validity and this speeds things up slightly
drawS3 <- function(x, setup, i) {
    ind <- getIndices(setup)[[i]]  # indices for i-th sample
    res <- x[ind, , drop=FALSE]
    prob <- getProb(setup)
    if(length(prob) > 0) res$.weight <- 1/(prob[ind])
    res
}
