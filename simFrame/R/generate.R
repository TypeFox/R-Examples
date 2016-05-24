# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

setMethod("generate",
    signature(control = "DataControl"),
    function(control) {
        # initializations
        size <- getSize(control)
        distribution <- getDistribution(control)
        dots <- getDots(control)
        nam <- getColnames(control)
        # generate data
        values <- do.call(distribution, c(size, dots))
        if(is.null(dim(values)) && is.null(nam)) nam <- "V1" 
        values <- as.data.frame(values)
        if(!is.null(nam)) {
            p <- ncol(values)
            if(length(nam) != p) {
                stop(gettextf("'names' must be a vector of length %i", p))
            }
            names(values) <- nam
        }
        values
    })

setMethod("generate", 
    signature(control = "character"), 
    function(control, ...) {
        if(length(control) != 1) {
            stop("'control' must specify exactly one ", 
                "class extending \"VirtualDataControl\"")
        }
        if(!extends(control, "VirtualDataControl")) {
            stop(gettextf("\"%s\" does not extend \"VirtualDataControl\"", 
                    control))
        }
        generate(new(control, ...))
    })

setMethod("generate",
    signature(control = "missing"),
    function(control, ...) {
        generate(DataControl(...))
    })
