###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# extractreplace.R
# Author: Francisco
###############################################################################
# DESCRIPTION: A set of utility functions for extracting and replacing data in COPViews and BLViews objects
# KEYWORDS: hplot
###############################################################################

viewMatrix <- function(views, dropZeroColumns = TRUE) {
    .assertClass(views, "BLViews")
    P <- views@P
    if(dropZeroColumns) {
        isZeroColumn <- apply(P==0, 2, all)
        P <- P[,!isZeroColumn, drop = FALSE]
    }
    cbind(P, "q" = views@qv)
}

"PMatrix<-" <- function(views, value)
{
    stopifnot(nrow(views@P) >= nrow(value)) 
#    stopifnot()
    views
}

PMatrix <- function(views)
{
   .assertClass(views, c("BLViews", "COPViews"))
    # TODO: rename "P" to "pick" in BLViews
    if(class(views) %in% "BLViews")
        views@P
    
    else
        views@pick
}

assetSet <- function(views)
{
    views@assets
}

"qv<-" <- function(views, value)
{
    .assertClass(views, "BLViews")
    if(length(value) != length(views@qv))
    {
        warning("Vector qv is of incorrect length, will not replace")
        return(views)   
    }
    views@qv <- value
    views
}

"confidences<-" <- function(views, value) {
    .assertClass(views, c("BLViews", "COPViews"))
    if(length(value) != length(views@confidences))
    {
        warning("value is of incorrect length, will not replace")
        return(views)   
    }
    views@confidences <- value
    views
}   


confidences <- function(views)
{
    .assertClass(views, c("BLViews", "COPViews"))
    views@confidences
}

# TODO: unit test

posteriorMeanCov <- function(posterior)
{
	.assertClass(posterior, "BLResult")
	list("covariance" = posterior@posteriorCovar, "mean" = posterior@posteriorMean)
	
}



## Extracts the matrix of posterior simulations from a COPPosterior object.
## Return A matrix with named columns.
posteriorSimulations <- function(posterior)
{
	.assertClass(posterior, "COPResult")
	posterior@posteriorSims
}

numSimulations <- function(posterior)
{
	.assertClass(posterior, "COPResult")
	nrow(posterior@posteriorSims)
}

priorViews <- function(posterior)
{
	.assertClass(posterior, "COPResult")
	posterior@views
}