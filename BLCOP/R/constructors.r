
###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# COPViews
# Author: Francisco
###############################################################################
# DESCRIPTION: Constructor function for the COPViews object
# KEYWORDS: utilities
###############################################################################

COPViews <- function
(
    pickMatrix,                   # View matrix
    viewDist,                     # list of marginal distributions of views
    confidences,                  # vector confidences in views
    assetNames                    # names of assets in one's "universe"
)
{
    if(is.null(colnames(pickMatrix)))
        colnames(pickMatrix) <- assetNames
    new("COPViews", pick = pickMatrix, viewDist = viewDist, 
        confidences = confidences, assets = assetNames)
}


###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# distribution
# Author: Francisco
###############################################################################
# DESCRIPTION: Constructor function for a distribution class object
# KEYWORDS: utilities
###############################################################################

distribution <- function
(
    RName,     # string holding the "R name" / suffix for an R distribution.  e.g. for a normal distribution this would be "norm"
    ...        # additional parameters for defining a distribution
)
{
    # construct the names of the associated distribution functions
    distFunctions <- paste(c("r" ,"d", "p", "q"), RName, sep = "")
    
    # at least the sampling function must exist for COP to be used
    if(!exists(distFunctions[1]))
        stop("Sampling function for this distribution does not exist!")
    
    
    samplingFun <- try(match.fun(distFunctions[1]))
    
    if("try-error" %in% class(samplingFun))
        stop(paste(samplingFun, "seems not to be a function!" ))
    if(any(!exists(distFunctions[-1])))
        warning("Some functions associated to this distribution could not be found")
    
    # extract the parameters passed in through the ellipsis as a named vector.  Then check that the sampling function actually accepts
    # all of these parameters 
    distParams <- as.numeric(list(...))
    names(distParams) <- names(list(...))
    expectedDistParams <- names(formals(distFunctions[1]))
    stopifnot(all(names(distParams) %in% expectedDistParams))
    
    new("distribution", RName = RName, parameters = distParams)
}

###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# mvdistribution
# Author: Francisco
###############################################################################
# DESCRIPTION: Constructor function for the multivariate distribution class object
# KEYWORDS: utilities
###############################################################################


mvdistribution <- function
(
    RName,          #string holding the "R name" / suffix for an R distribution.  e.g. for a normal distribution this would be "norm"
    ...             # additional parameters for defining a distribution
)
{
    distFunctions <- paste(c("r" ,"d", "p", "q"), RName, sep = "")
    if(!exists(distFunctions[1]))
        stop("Sampling function for this distribution does not exist!")
    if(any(!exists(distFunctions[-1])))
        warning("Some functions associated to this distribution could not be found")
    
    samplingFun <- try(match.fun(distFunctions[1]))    
    if("try-error" %in% class(samplingFun))
        stop(paste(samplingFun, "seems not to be a function!" ))
    
    distParams <- list(...)
    names(distParams) <- names(distParams)
    expectedDistParams <- names(formals(distFunctions[1]))
    stopifnot(all(names(distParams) %in% expectedDistParams))
    new("mvdistribution", RName = RName, parameters = distParams)
}

###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# BLViews
# Author: Francisco
###############################################################################
# DESCRIPTION: Constructor function for a BLViews object.
# KEYWORDS: datagen
###############################################################################
BLViews <- function
(                          
    P,                     # Pick matrix 
    q,                     # vector of "q" values in the Black-Litterman model
    confidences,           # vector of confidences in views
    assetNames             # names of assets
)
{
    
    if(is.null(colnames(P)))
        dimnames(P) <- list(rownames(P), assetNames)
    names(q) <- NULL
    new("BLViews", P = P, qv = q, confidences = confidences, assets = assetNames)
}
