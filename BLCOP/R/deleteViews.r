###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# deleteViews
# Author: Francisco
###############################################################################
# DESCRIPTION: Generic function that deletes a vector of views (which correspond to rows of the pick matrix) from a view object.  returns
# the view object with the views deleted
# KEYWORDS: manip, utilities
###############################################################################


deleteViews <- function(views, viewsToDel) 
{ 
    invisible(NULL)
}

setGeneric("deleteViews")

deleteViews.BLViews <- function(views, viewsToDel) 
{
    if(any(viewsToDel > nrow(views@P)))
    {    
        warning("Attempting to delete some non-existent views, ignoring")
        viewsToDel <- viewsToDel[viewsToDel < nrow(views@P)]
    }
    views@P <- views@P[-viewsToDel,,drop=FALSE]
    views@qv <- views@qv[-viewsToDel]
    views@confidences <- views@confidences[-viewsToDel]
    views
}

setMethod("deleteViews", signature(views = "BLViews"), deleteViews.BLViews)

deleteViews.COPViews <- function(views, viewsToDel)
{
    if(any(viewsToDel > nrow(views@pick)))
    {    
        warning("Attempting to delete some non-existent views, ignoring")
        viewsToDel <- viewsToDel[viewsToDel < nrow(views@pick)]
    }
    views@pick <- views@pick[-viewsToDel,,drop=FALSE]
#viewDist = "list"
    views@viewDist <- views@viewDist[-viewsToDel]
    views@confidences <- views@confidences[-viewsToDel]
    views
}
setMethod("deleteViews", signature(views = "COPViews"), deleteViews.COPViews)