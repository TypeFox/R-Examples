###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# createBLViews
# Author: Francisco
###############################################################################
# DESCRIPTION: Creates a set of B-L views with a GUI interface.  Currently the only interface that can be used is the data editor
# KEYWORDS: datagen
###############################################################################

createBLViews <- function
(
  allAssets,                   # A vector of strings holding the names of all of the assets in one's universe
  numAssetViews = 1, #         # Number of views to form
  assetSubset = NULL,          # subset of one's universe that one actually wants to form views on
  mode = c("editor", "Window") # GUI to use, only editor used at the moment
)
{
  mode <- match.arg(mode)
  .createBLViews.Editor(allAssets, numAssetViews, assetSubset)
}

###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# createViews
# Author: Francisco
###############################################################################
# DESCRIPTION: Creates a set of B-L views using the data editor as a GUI
# KEYWORDS: datagen
###############################################################################


.createBLViews.Editor <- function 
(
  allAssets,
  numViews = 1,
  stockSubset = NULL
)
{
  #extract 
  if( is.null(stockSubset) )                                                                                        
    stockSubset <- allAssets
  if(length(stockSubset) < numViews)
    stop("The number of views to be formed should be greater than or equal to the number of stocks")
  viewsMatInit <- matrix(0, nrow = numViews, ncol = length(stockSubset) + 2, 
                  dimnames = list(NULL, c( stockSubset, "q", "confidence")))
  viewsMat <- edit(viewsMatInit)
  
  P <- viewsMat[, -((ncol(viewsMat)-1):ncol(viewsMat)), drop = FALSE] 
  qv <- viewsMat[,ncol(viewsMat) - 1 ]
  conf <- viewsMat[, ncol(viewsMat)]
  
  views <- try(BLViews(P,q = qv, conf, assetNames = allAssets))
  if(inherits(views, "try-error"))
    stop("Incorrectly entered views, unable to initialize object")
  views
}

updateBLViews <- function (
  views,
  includeNullViews = FALSE,
  numNewViews = 0,
   assets = NULL
)
{
    
    if(!is.null(assets))  
        assets <- assetSet(views)
      
  # extract the "P" matrix and then remove those stocks for which we have not formed
  # any views at all, i.e. those with only 0 column entries      
    if(!includeNullViews)
        P <- .removeZeroColumns(views@P)
    else
        P <- views@P
    temp <- colnames(P)
    
    newAssets <- setdiff(assets, temp)
    addCols <- matrix(0, nrow = nrow(P), ncol = length(newAssets), dimnames = list(NULL, newAssets))
     
  # aggregate pre-existing views    
    P <- cbind(P, addCols)
    if(numNewViews > 0)
        P <- rbind(P, matrix(0, ncol = ncol(P), nrow = numNewViews))
           
    qv <- c(views@qv, rep(0, numNewViews))
    confidences <- c(views@confidences, rep(0, numNewViews))
  
    viewsMatInit <- cbind(P, "qv" = qv, "confidences" = confidences)  
    viewsMat <- edit(viewsMatInit)
   
    v <- viewsMat[, -ncol(viewsMat)] 
    conf <- viewsMat[, ncol(viewsMat)]
  
    BLViews(P = v[, -ncol(v)], q = v[,ncol(v)],  confidences = conf, assetNames = colnames(P))
}