###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# createCOPViews
# Author: Francisco
###############################################################################
# DESCRIPTION: Creates _part_ of a set of COP views with a GUI interface.  Currently the only interface that can be used is the 
# data editor, and no method exists for specifying the views graphically.  Only the pick matrices and confidences can be set with this.
# KEYWORDS: datagen
###############################################################################


createCOPViews <- function(
  allAssets,                       # Asset universe
  numAssetViews = 1, #             # Number of views to form
  assetSubset = NULL,              # Subset of asset universe to form views on
  mode = c("editor", "Window")     # Currently unused.
)
{
#  stopifnot((is.null(stockSet) && numFactorViews == 0) || ( !is.null(stockSet) && numFactorViews > 0))
  mode <- match.arg(mode)
  .createCOPViews.Editor(allAssets, numAssetViews, assetSubset)
}


.createCOPViews.Editor <- function (
  allAssets,
  numViews = 1,
  assetSubset = NULL
)
{
  DEFAULTCONFIDENCE <- 1/10000
  #extract 
  if( is.null(assetSubset) )                                                                                        
    assetSubset <- allAssets
  if(length(assetSubset) < numViews)
    stop("The number of views to be formed should be greater than or equal to the number of assets")
  viewsMatInit <- matrix(0, nrow = numViews, ncol = length(assetSubset) + 1, 
                  dimnames = list(NULL, c( assetSubset, "confidence")))
  viewsMat <- edit(viewsMatInit)
  
  # the pick matrix excludes the last column, which has the confidences
  P <- viewsMat[, -ncol(viewsMat), drop = FALSE] 
  conf <- viewsMat[, "confidence"]
  viewDist <- lapply(rep("norm", numViews), distribution, mean = 0, 
  sd = DEFAULTCONFIDENCE)

  views <- try(COPViews(P, viewDist, conf,assetNames = allAssets))
  if(inherits(views, "try-error"))
    stop("Incorrectly entered views, unable to initialize object")
  views
}