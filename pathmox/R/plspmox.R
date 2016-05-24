#' @title PLS Path Model of a node from a PATHMOX or TECHMOX tree
#' 
#' @description
#' Calculates a PLS Path Model on a selected node from a pathmox or techmox tree
#' 
#' @details
#' Performs a PLS-PM analysis with the elements contained in \code{node}, by
#' calling the function \code{\link{plspm}}. The rest of the parameters to
#' perform the PLS-PM analysis (i.e. \code{inner}, \code{outer}, \code{modes},
#' \code{scheme}, \code{scaled}, \code{tol}, \code{iter}) are inherited from the
#' object in argument \code{pls}. \cr When the object \code{pls} does not
#' contain a data matrix (i.e. \code{pls$data=NULL}), the user must provide the
#' data matrix or data frame in \code{X}.
#' 
#' @param pls An object of class \code{"plspm"}
#' @param treemox An object of class \code{"treemox"}
#' @param X Optional argument for data table
#' @param node An integer value indicating the number of the node. Must be an
#' integer larger than 1.
#' @param boot.val A logical value indicating whether bootstrap validation is
#' performed (\code{FALSE} by default).
#' @param br An integer indicating the number bootstrap resamples. Used only
#' when \code{boot.val=TRUE}.
#' @param dataset A logical value indicating whether the data matrix should be
#' included in the list of results (\code{FALSE} by default).
#' @return An object of class \code{"plspm"}
#' @author Gaston Sanchez
#' @seealso See Also as \code{\link{plspm}}, \code{\link{pathmox}}
#' @export
#' @examples
#' 
#'  \dontrun{
#'  ## example of PLS-PM in customer satisfaction analysis
#'  ## model with seven LVs and reflective indicators
#'  data(csimobile)
#'  
#'  # select manifest variables
#'  data_mobile = csimobile[,8:33]
#'  
#'  # define path matrix (inner model)
#'  IMAG = c(0, 0, 0, 0, 0, 0, 0)
#'  EXPE = c(1, 0, 0, 0, 0, 0, 0)
#'  QUAL = c(0, 1, 0, 0, 0, 0, 0)
#'  VAL = c(0, 1, 1, 0, 0, 0, 0)
#'  SAT = c(1, 1, 1, 1, 0, 0, 0)
#'  COM = c(0, 0, 0, 0, 1, 0, 0)
#'  LOY = c(1, 0, 0, 0, 1, 1, 0)
#'  mob_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, COM, LOY)
#'  
#'  # blocks of indicators (outer model)
#'  mob_blocks = list(1:5, 6:9, 10:15, 16:18, 19:21, 22:24, 25:26)
#'  mob_modes = rep("A", 7)
#'  
#'  # apply plspm
#'  mob_pls = plspm(data_mobile, mob_path, mob_blocks, 
#'                  modes = mob_modes, scheme="factor", scaled=FALSE)
#'
#'  # re-ordering those segmentation variables with ordinal scale (Age and Education)
#'  csimobile$Education = factor(csimobile$Education, 
#'      levels=c("basic","highschool","university"),
#'      ordered=TRUE)
#'  
#'  # select the segmentation variables
#'  seg_vars = csimobile[,1:7]
#'  
#'  # Pathmox Analysis
#'  mob_pathmox = pathmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)
#'
#'  # get PLS-PM of nodes 2 and 3
#'  node2 = plspmox(mob_pls, mob_pathmox, node=2)
#'  node3 = plspmox(mob_pls, mob_pathmox, node=3)
#'  }
#'
plspmox <-
function(pls, treemox, X=NULL, node, boot.val=FALSE, br=NULL, dataset=FALSE)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm")
    stop("\nAn object of class 'plspm' was expected")
  if (class(treemox) != "treemox")
    stop("\nArgument 'treemox' must be an object of class 'treemox'")
  if (nrow(pls$scores) != treemox$MOX$Size[1]) 
    stop("\n'pls' and 'treemox' have different number of observations")
  if (!is.null(X)) # if X available
  {
    if (is.null(pls$data))
    {
      if (is_not_tabular(X))
        stop("\nparameter 'X' must be a numeric matrix or data frame.")
      if (nrow(X) != nrow(pls$scores))
        stop("\n'pls' and 'X' have different number of rows.")
    }
  } else { # if no X
    if (is.null(pls$data)) 
      stop("\n'X' is missing. No dataset available.")
  }
  if (missing(node))
    stop("\nargument 'node' (number of node) is missing")
  if (mode(node)!="numeric" || length(node)!=1 || node<=1 || (node%%1)!=0)
    stop("\nInvalid number of 'node'. Must be an integer larger than 1")
  if (length(which(treemox$MOX$Node==node)) == 0)
    stop("\nInvalid number of 'node'")
  
  # =======================================================
  # inputs setting
  # =======================================================  
  IDM <- pls$model$IDM
  blocks <- pls$model$blocks
  modes <- pls$model$specs$modes
  scheme <- pls$model$specs$scheme
  scaled <- pls$model$specs$scaled
  tol <- pls$model$specs$tol
  maxiter <- pls$model$specs$maxiter

  # data matrix DT
  if (!is.null(pls$data)) {
    DT <- pls$data
  } else {
    # building data matrix 'DT' from 'X'
    DT = get_manifests(X, blocks)
  }
  # get node and its observations
  list.nodes <- treemox$list.nodes
  node <- which(treemox$MOX$Node == node) - 1
  node.obs <- list.nodes[[node]]
  DT.node <- DT[node.obs,]
  
  # apply plspm
  res <- plspm(DT.node, IDM, blocks, scaling = NULL, modes, 
               scheme, scaled, boot.val=boot.val, 
               br=br, tol=tol, maxiter=maxiter, dataset=FALSE)
  return(res)
}

