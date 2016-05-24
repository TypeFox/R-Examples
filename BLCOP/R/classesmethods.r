######################
# validity functions
######################

.BLViews.valid <- function(object)
{
    numViews <- c(nrow(object@P), length(object@qv), length(object@confidences))
    if(length(unique(numViews)) != 1)
        return("Inconsistent number of views implied\n")
    if(any(object@confidences < 0))
        return("Negative confidences")
    if(!setequal(colnames(object@P), object@assets))
        return("asset names not consistent with P's column names")
    return(TRUE)
}

.BLResult.valid <- function(object)
{
    if(length(object@posteriorMean) != length(object@priorMean))
    {
        return(FALSE)
    }
    if(!all(dim(object@priorCovar) == dim(object@posteriorMean)))
        return(FALSE)
    TRUE
}

.COPViews.valid <- function(object)
{
    # All of these quantities should be equal
    numViews <- c(nrow(object@pick), length(object@viewDist), 
    length(object@confidences))
    #check if any don't match
    if(length(unique(numViews)) != 1)
        return("Inconsistent number of views implied!\n ")
    
    if(any(object@confidences < 0 | object@confidences > 1))
        return("Confidences must lie between 0 and 1 \n")
    
    if(!setequal(colnames(object@pick), object@assets))
        return("asset names not consistent with pick matrix's column names \n")
    return(TRUE)
}

.COPResult.valid <- function(object) 
{    
    if(length(object@views@assets) != ncol(object@posteriorSims) )
        return(FALSE)
    return(TRUE)
}

setClass("BLViews", representation(P = "matrix", qv = "numeric", confidences = "numeric", 
    assets = "character"), validity = .BLViews.valid)
setClass("BLResult", representation(views = "BLViews", tau = "numeric", priorMean = "numeric", priorCovar = "matrix", 
    posteriorMean = "numeric", posteriorCovar = "matrix", kappa = "numeric"), validity = .BLResult.valid)
setClass("distribution", representation(RName = "character", parameters = "numeric"))
setClass("mvdistribution", representation(RName = "character", parameters = "list"))
setClass("COPViews", representation(pick = "matrix", viewDist = "list",confidences = "numeric", assets = "character" ),
    validity = .COPViews.valid )
setClass("COPResult", representation(views = "COPViews", marketDist = "mvdistribution", posteriorSims = "matrix"),
            validity = .COPResult.valid)