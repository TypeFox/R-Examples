# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## convenience wrapper for 'setup'

simSample <- function(x, design = character(), grouping = character(), 
        collect = FALSE, fun = srs, size = NULL, 
        prob = NULL, ..., k = 1) {
    # define control object
    control <- SampleControl(design=design, grouping=grouping, 
        collect=collect, fun=fun, size=size, prob=prob, dots=list(...), k=k)
    # call 'setup'
    res <- setup(x, control)
    call <- match.call()
    setCall(res, call)
    res
}
