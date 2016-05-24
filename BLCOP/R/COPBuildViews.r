
###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# addCOPViews
# Author: Francisco
###############################################################################
# DESCRIPTION: Adds views specified in a pick matrix, set of confidences and view distributions to a pre-existing
# view object
# KEYWORDS: datagen, manip
###############################################################################

addCOPViews <- function
(
    pickMatrix,                           #  pick matrix
    viewDist,                             #  list of distribution of views
    confidences,                          #  vector of confidences in views
    views                                 #  pre-existing object to add views to
)
{
 
    
    sNames <- colnames(pickMatrix)
    if(is.null(sNames))
        stop("Missing asset names in pick matrix, cannot input views")
    
    assets <- assetSet(views)
    positions <- match(sNames, assets)
    if(any(is.na(positions)))
        stop("Some asset names in the pick matrix did not have matches to list of provided assets")
    
    # construct a new pick matrix to incorporate old and new views
    P <- matrix(0, ncol = length(assets), nrow = nrow(pickMatrix), dimnames = list(NULL, assets))
    # insert the new views
    P[, positions ] <- pickMatrix
 
    new("COPViews", "pick" = rbind(views@pick, P), viewDist = c(views@viewDist, viewDist),
                  "confidences" = c(views@confidences, confidences), "assets" = assets)    
}