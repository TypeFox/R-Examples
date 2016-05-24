
# $Id: Classes.R 532 2014-02-01 08:56:55Z thothorn $

### Linear statistic with expectation and covariance
setClass(Class = "LinStatExpectCovar",
    representation = representation(
        linearstatistic = "numeric",
        expcovinf = "ExpectCovarInfluence"
    ),
    contains = "ExpectCovar"
)

### Memory for C_svd
setClass(Class = "svd_mem",
    representation = representation(
        method = "character",
        jobu   = "character",
        jobv   = "character",
        u      = "matrix",
        v      = "matrix",
        s      = "numeric",
        p      = "integer"
    )
)

### with Moore-Penrose inverse of the covariance matrix
setClass(Class = "LinStatExpectCovarMPinv",
    representation = representation(
        MPinv  = "matrix",   
        rank   = "numeric",
        svdmem = "svd_mem"
    ), 
    contains = "LinStatExpectCovar"
)

################ Memory Classes #####################

setClass(Class = "TreeFitMemory",
    representation = representation(
        expcovinf         = "ExpectCovarInfluence",
        expcovinfss       = "ExpectCovarInfluence",
        linexpcov2sample  = "LinStatExpectCovar",
        weights           = "numeric",
        varmemory         = "list",
        dontuse           = "logical",
        dontusetmp        = "logical",
        splitstatistics   = "numeric"
    ), 
    validity = function(object) {
        ni <- length(dontuse)
        length(varmemory) == ni && length(dontusetmp) == ni
    }
)


##############  Tree Classes  ######################

setClassUnion("df_OR_list", c("data.frame", "list"))

setClass(Class = "VariableControl",
    representation = representation(
        teststat = "factor",
        pvalue   = "logical",
        tol      = "numeric",
        maxpts   = "integer",
        abseps   = "numeric",
        releps   = "numeric"
    ),
    prototype = list(
        teststat = factor("max", levels = c("max", "quad")),
        pvalue   = as.logical(TRUE),
        tol      = as.double(1e-10),
        maxpts   = as.integer(25000),
        abseps   = as.double(1e-4),
        releps   = as.double(0.0)
    )
)

setClass(Class = "SplitControl",
    representation = representation(
        minprob      = "numeric",
        minsplit     = "numeric",
        minbucket    = "numeric",
        tol          = "numeric",
        maxsurrogate = "integer"
    ),
    prototype = list(minprob = as.double(0.01), 
                     minsplit = as.double(20), 
                     minbucket = as.double(7), 
                     tol = as.double(1e-10), 
                     maxsurrogate = as.integer(0)
    ),
    validity = function(object) {
        if (any(c(object@minsplit, object@minbucket, 
                  object@tol, object@maxsurrogate) < 0)) {
            warning("no negative values allowed in objects of class ", 
                    sQuote("SplitControl"))
            return(FALSE)
        }
        if (object@minprob < 0.01 || object@minprob > 0.99) {
            warning(sQuote("minprob"), " must be in (0.01, 0.99)")
            return(FALSE)
        }
        return(TRUE)
    }
)

setClass(Class = "GlobalTestControl",
    representation = representation(
        testtype     = "factor",
        nresample    = "integer",
        randomsplits = "logical",
        mtry         = "integer",
        mincriterion = "numeric"
    ),
    prototype = list(
        testtype = factor("Bonferroni", 
            levels = c("Bonferroni", "MonteCarlo", "Aggregated", 
                       "Univariate", "Teststatistic")),
        nresample = as.integer(9999),
        randomsplits = as.logical(FALSE),
        mtry = as.integer(0),
        mincriterion = as.double(0.95)
    ),
    validity = function(object) {
        if (object@mincriterion < 0) {
            warning(sQuote("mincriterion"), " must not be negative")
            return(FALSE)
        }
        if (any(object@mtry < 0)) {
            warning(sQuote("mtry"), " must be positive")
            return(FALSE)
        }
        if (object@nresample < 100) {
            warning(sQuote("nresample"), " must be larger than 100")
            return(FALSE)
        }
        return(TRUE)
    },
)

setClass(Class = "TreeGrowControl",
    representation = representation(
        stump           = "logical",
        maxdepth        = "integer",
        savesplitstats  = "logical"
    ),
    prototype = list(stump = as.logical(FALSE), 
                     maxdepth = as.integer(0), 
                     savesplitstats = as.logical(TRUE)),
    validity = function(object) {
        if (object@maxdepth < 0) {
            warning(sQuote("maxdepth"), " must be positive")
            return(FALSE)
        }
        return(TRUE)
    }
)

setClass(Class = "TreeControl",
    representation = representation(
        varctrl   = "VariableControl",
        splitctrl = "SplitControl",
        gtctrl    = "GlobalTestControl",
        tgctrl    = "TreeGrowControl"
    ),
    prototype = list(varctrl = new("VariableControl"),
                     splitctrl = new("SplitControl"),
                     gtctrl = new("GlobalTestControl"),
                     tgctrl = new("TreeGrowControl")
    ),
    validity = function(object) {
        (validObject(object@varctrl) && 
        validObject(object@splitctrl)) &&
        (validObject(object@gtctrl) &&
        validObject(object@tgctrl))
    }
)

setClass(Class = "ForestControl",
    representation = representation(
        ntree    = "integer",
        replace  = "logical",
        fraction = "numeric",
        trace    = "logical"),
    contains = "TreeControl",
    validity = function(object) {
        if (object@ntree < 1) {
            warning(sQuote("ntree"), " must be equal or greater 1")
            return(FALSE)
        }
        if (object@fraction < 0.01 || object@fraction > 0.99) {
            warning(sQuote("fraction"), " must be in (0.01, 0.99)")
            return(FALSE)
        }
        return(TRUE)
    }
)

setClass(Class = "VariableFrame",
    representation = representation(
        variables       = "df_OR_list", 
        transformations = "list", 
        is_nominal      = "logical", 
        is_ordinal      = "logical",
        is_censored     = "logical",
        ordering        = "list", 
        levels          = "list", 
        scores          = "list",
        has_missings    = "logical", 
        whichNA         = "list",
        nobs            = "integer",
        ninputs         = "integer")
)

setClass(Class = "ResponseFrame",
    representation = representation(
        test_trafo = "matrix",
        predict_trafo = "matrix"
    ), contains = "VariableFrame"
)   

setClass(Class = "LearningSample",
    representation = representation(
        responses = "ResponseFrame",
        inputs    = "VariableFrame",
        weights   = "numeric",
        nobs      = "integer",
        ninputs   = "integer"
    )
)

setClass(Class = "LearningSampleFormula",
    representation = representation(
        menv      = "ModelEnv"
    ), contains = "LearningSample"
)

### the tree structure itself is a list, 
### and we need to make sure that the tree slot excepts
### the S3 classes. 
setClass(Class = "SplittingNode", contains = "list")
setClass(Class = "TerminalNode", contains = "list")
setClass(Class = "TerminalModelNode", contains = "list")
setClass(Class = "orderedSplit", contains = "list")
setClass(Class = "nominalSplit", contains = "list")

### and we don't want to see warnings that class `Surv'
### (S3 method in `survival') is unknown
setClass(Class = "Surv", contains = "list")


### A class for partitions induced by recursive binary splits
setClass(Class = "BinaryTreePartition",
    representation = representation(
        tree     = "list",          # the basic tree structure as (named or
                                    # unnamed) list
        where    = "integer",       # the nodeID of the observations in the
                                    # learning sample
        weights  = "numeric"         # the weights in the root node
    ),
)

### A class for binary trees   
setClass(Class = "BinaryTree", 
    representation = representation(
        data                = "ModelEnv",
        responses           = "VariableFrame", # a list of response `variables'
                                               # for computing predictions
        cond_distr_response = "function",      # predict distribtion
        predict_response    = "function",      # predict responses
        prediction_weights  = "function",      # prediction weights
        get_where           = "function",      # node numbers
        update              = "function"       # update weights
    ),
    contains = "BinaryTreePartition"
)

### A class for random forest  
setClass(Class = "RandomForest", 
    representation = representation(
        ensemble            = "list",
        where               = "list",
        weights             = "list",
        initweights         = "numeric",
        data                = "ModelEnv",
        responses           = "VariableFrame", # a list of response `variables'
                                               # for computing predictions
        cond_distr_response = "function",      # predict distribtion
        predict_response    = "function",      # predict responses
        prediction_weights  = "function",      # prediction weights
        get_where           = "function",      # node numbers
	update              = "function"       # update weights
    )
)

