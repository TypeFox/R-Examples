### S4 class definitions for lsmeans package


### ref.grid object -- for a reference grid
setClass("ref.grid", representation (
    model.info = "list",
    roles = "list",
    grid = "data.frame", 
    levels = "list",
    matlevs = "list",
    linfct = "matrix",
    bhat = "numeric",
    nbasis = "matrix",
    V = "matrix",
    dffun = "function",
    dfargs = "list",
    misc = "list",
    post.beta = "matrix"
))
# Note: misc will hold various extra params,
# including at least the following req'd by the summary method
#   estName: column name for the estimate in the summary ["prediction"]
#   infer: booleans (CIs?, tests?)  [(FALSE,FALSE)]
#   level: default conf level [.95]
#   adjust: default adjust method ["none"]
#   famSize: number of means in family

### lsmobj class -- almost trivial ext of ref.grid, structurally
# But origin can be very different from those of a reference grid
# In general its 'grid' will correspond to some set of 
# linear functions of grid points
setClass("lsmobj", contains="ref.grid")

