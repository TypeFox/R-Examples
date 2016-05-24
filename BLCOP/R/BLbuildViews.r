
###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# addBLViews
# Author: Francisco
###############################################################################
# DESCRIPTION: Adds a new set of views
# KEYWORDS: datagen
# TODO: Add the ability to add entirely new assets
###############################################################################

addBLViews <- function
(
    pickMatrix,  # pickMatrix matrix.
    q, 	 # mean vector
    confidences,   # vector of confidences in the new views
    views         # pre-existing views to add to
)
{

    # try to force the view input to be a matrix just in case if its a vector
    if(class(pickMatrix) == "numeric")
        pickMatrix <- matrix(pickMatrix, nrow = 1, ncol = length(pickMatrix), dimnames = list(NULL, names(pickMatrix)) )  
        
    
    if(is.null(colnames(pickMatrix)))
    {    
        warning("Missing asset names in the pickMatrix matrix, assigning them automatically")
        dimnames(pickMatrix) <- list(rownames(pickMatrix), assetSet(views)[1:ncol(pickMatrix)])
    }
    sNames <- colnames(pickMatrix)
    assetNames <- assetSet(views)
    
    # find the indices of the names of the assets which occur in the new views
    # within the vector of asset names of the already-existing views
    positions <- match(sNames, assetNames)
    if(any(is.na(positions)))
        stop("Some asset names in the new views matrix did not have matches to the assetnames in the first object")
    
    P <- matrix(0, ncol = length(assetNames), nrow = nrow(pickMatrix), dimnames = list(NULL, assetNames))
    P[, positions ] <- pickMatrix

  # create a new set of views
    q <- c(views@qv, q)
    names(q) <- NULL
    views <- new("BLViews", "P" = rbind(views@P, P), "qv" = q,
                "confidences" = c(views@confidences, confidences), "assets" = assetNames)

}