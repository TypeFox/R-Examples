# S4 class definition for exemplar-based clustering
setClass("ExClust",
    representation = representation
    (
        l         = "numeric",
        sel       = "numeric",
        exemplars = "numeric",
        clusters  = "list",
        idx       = "numeric",
        sim       = "mMatrix",
        call      = "character"
    ),
    prototype = prototype
    (
        l         = 0,
        sel       = numeric(0),
        exemplars = numeric(0),
        clusters  = list(),
        idx       = numeric(0),
        sim       = matrix(nrow=0, ncol=0),
        call      = character(0)
    )
)

# S4 class definition for the result object of affinity propagation clustering
setClass("APResult",
    representation = representation
    (
        sweeps    = "numeric",
        it        = "numeric",
        p         = "numeric",
        netsim    = "numeric",
        dpsim     = "numeric",
        expref    = "numeric",
        netsimLev = "numeric",
        netsimAll = "numeric",
        dpsimAll  = "numeric",
        exprefAll = "numeric",
        idxAll    = "matrix"
    ),
    prototype = prototype
    (
        sweeps    = 0,
        it        = 0,
        p         = 0,
        netsim    = NaN,
        dpsim     = NaN,
        expref    = NaN,
        netsimLev = numeric(0),
        netsimAll = NaN,
        dpsimAll  = NaN,
        exprefAll = NaN,
        idxAll    = matrix(nrow=0, ncol=0)
    ),
    contains = "ExClust"
)

# S4 class definition for the result object of the aggExCluster algorithm
setClass("AggExResult",
    representation = representation
    (
        l             = "numeric",
        sel           = "numeric",
        maxNoClusters = "numeric",
        clusters      = "list",
        exemplars     = "list",
        merge         = "matrix",
        height        = "numeric",
        order         = "numeric",
        labels        = "character",
        sim           = "matrix",
        call          = "character"
    ),
    prototype = prototype
    (
        l             = 0,
        sel           = numeric(0),
        maxNoClusters = 0,
        clusters      = list(),
        exemplars     = list(),
        merge         = matrix(NA, 1, 1),
        height        = numeric(0),
        order         = numeric(0),
        labels        = c(),
        sim           = matrix(NA, 1, 1),
        call          = character(0)
    )
)
