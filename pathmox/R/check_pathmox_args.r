#' @title Check arguments for \code{pathmox} and \code{pathmox.fit}
#' 
#' @details
#' Internal function. \code{check_pathmox_args} is called by \code{pathmox}.
#' 
#' @param pls An object of class \code{"plspm"} returned by \code{\link{plspm}}.
#' @param EXEV A data frame of factors contaning the segmentation variables.
#' @param X Optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{pls}.
#' @param signif A numeric value indicating the significance threshold of the
#' F-statistic. Must be a decimal number between 0 and 1.
#' @param size numeric value indicating the minimum size of elements inside a
#' node.
#' @param deep integer indicating the depth level of the tree. Must be an
#' integer greater than 1.
#' @param tree whether the tree should be displayed
#' (\code{TRUE} by default).
#' @return list of validated arguments
#' @keywords internal
#' @export
check_pathmox_args <-
function(pls, EXEV, X, signif, size, deep, tree)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm")
    stop("\n'pls' must be an object of class 'plspm'")
  
  # check compatibility of 'pls$data' or 'X'
  if (!is.null(X)) # if X available
  {
    if (is.null(pls$data))
    {
      if (is_not_tabular(X))
        stop("\n'X' must be a numeric matrix or data frame")
      if (nrow(X) != nrow(pls$scores))
        stop("\n'pls' and 'X' have different number of rows")
      if (nrow(X) != nrow(EXEV))
        stop("\n'X' and 'EXEV' have different number of rows")
    }
  } else { 
    # if no X
    if (is.null(pls$data)) {
      stop("\n'X' is missing. No dataset available.")
    } else {
      if (nrow(pls$data) != nrow(EXEV)) 
        stop("\n'pls' and 'EXEV' have different number of rows")
    }
  }
  
  # check EXEV
  if (!is.data.frame(EXEV)) 
    stop("\n'EXEV' must be a data frame containing factors")
  for (j in 1:ncol(EXEV)) {
    if (!is.factor(EXEV[,j]))
      stop("\nOne or more columns in 'EXEV' are not factors")        
  }
  
  # check signif
  if (mode(signif)!="numeric" || length(signif)!=1 || signif<=0 || signif>=1)
  {
    cat("NOTICE: Invalid argument 'signif'. Default value 0.05 was used", "\n")
    signif <- 0.05
  }
  
  # check size
  if (mode(size)!="numeric" || length(size)!=1 || size<=0)
  {
    cat("NOTICE: Invalid argument 'size'. Default value 0.10 was used", "\n")
    size <- 0.10
  }        
  
  # check deep
  if (mode(deep)!="numeric" || length(deep)!=1 || deep<=1 || (deep%%1)!=0)
  {
    cat("NOTICE: Invalid argument 'deep'. Default value 2 was used", "\n")
    deep <- 2
  } 
  
  # check tree
  if (!is.logical(tree)) tree = TRUE
  # list with verified arguments
  list(signif = signif,
       size = signif,
       deep = deep,
       tree = tree)
}
