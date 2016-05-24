# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

setMethod("stratify", 
    signature(x = "data.frame", design = "BasicVector"), 
    function(x, design) {
        if(getSelectionLength(design) == 0) {
            stop("'design' must be specified")
        }
        tmp <- x[, design, drop=FALSE]
        tab <- as.data.frame(table(tmp))  # table of frequencies
        size <- tab[, ncol(tab)]  # strata sizes
        tab <- tab[, -ncol(tab), drop=FALSE]  # strata legend
        names(tab) <- names(tmp)
        spl <- split(1:nrow(tmp), tmp)  # split indices according to design
        names(spl) <- NULL
        nr <- 1:length(size)  # strata numbers
        values <- unsplit(nr, tmp)  # now we know stratum of each observation
        Strata(values=values, split=spl, 
            design=names(tmp), nr=nr, legend=tab, size=size)
    })


## utilities

# get information about strata as data.frame
setMethod("getStrataLegend",
    signature(x = "data.frame", design = "BasicVector"),
    function(x, design) {
        tab <- getStrataTable(x, design)
        tab[, -ncol(tab), drop=FALSE]
    })

# get list of indices in the strata
setMethod("getStrataSplit",
    signature(x = "data.frame", design = "BasicVector"),
    function(x, design, USE.NAMES = TRUE) {
        res <- split(1:nrow(x), x[, design])
        if(!USE.NAMES) names(res) <- NULL
        res
    })

# get contincency table as data.frame
setMethod("getStrataTable",
    signature(x = "data.frame", design = "BasicVector"),
    function(x, design) {
        tmp <- x[, design, drop=FALSE]
        ans <- as.data.frame(table(tmp))
        names(ans) <- c(names(tmp), "Size")
        ans
    })

# get strata sizes
setMethod("getStratumSizes",
    signature(x = "data.frame", design = "BasicVector"),
    function(x, design, USE.NAMES = TRUE) {
        spl <- getStrataSplit(x, design)
        getStratumSizes(spl, USE.NAMES=USE.NAMES)
    })

setMethod("getStratumSizes",
    signature(x = "list", design = "missing"),
    function(x, design, USE.NAMES = TRUE) {
        sapply(x, length, USE.NAMES=USE.NAMES)
    })

# get number of stratum for each observation
setMethod("getStratumValues",
    signature(x = "data.frame", design = "BasicVector", split = "missing"),
    function(x, design, split) {
        spl <- getStrataSplit(x, design)
        getStratumValues(x, design, spl)
    })

setMethod("getStratumValues",
    signature(x = "data.frame", design = "BasicVector", split = "list"),
    function(x, design, split) {
        unsplit(1:length(split), x[, design])
    })
