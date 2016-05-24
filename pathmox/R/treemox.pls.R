#' @title PLS-PM results of terminal nodes from a PATHMOX or TECHMOX tree
#' 
#' @description
#' Calculates basic PLS-PM results for the terminal nodes of PATHMOX and TECHMOX
#' trees
#' 
#' @details
#' The argument \code{pls} must be the same used for calculating the
#' \code{treemox} object.  When the object \code{pls} does not contain a data
#' matrix (i.e. \code{pls$data=NULL}), the user must provide the data matrix or
#' data frame in \code{X}.
#' 
#' @param pls An object of class \code{"plspm"} returned by \code{\link{plspm}}.
#' @param treemox An object of class \code{"treemox"} returned by
#' \code{\link{pathmox}} or \code{\link{techmox}}.
#' @param X Optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{pls}.
#' @return An object of class \code{"treemox.pls"}. Basically a list with the
#' following results:
#' @return \item{weights}{Matrix of outer weights for each terminal node}
#' @return \item{loadings}{Matrix of loadings for each terminal node}
#' @return \item{paths}{Matrix of path coefficients for each terminal node}
#' @return \item{r2}{Matrix of r-squared coefficients for each terminal node}
#' @author Gaston Sanchez
#' @seealso \code{\link{pathmox}}, \code{\link{techmox}},
#' \code{\link{plot.treemox}}.
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
#'  mob_pls = plspm(data_mobile, mob_path, mob_blocks, modes = mob_modes, 
#'                  scheme = "factor", scaled = FALSE)
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
#'  # applying function treemox.pls
#'  mob_nodes = treemox.pls(mob_pls, mob_pathmox)
#'  
#'  # comparative barplots
#'  plot(mob_nodes)
#'  }
#'
treemox.pls <-
function(pls, treemox, X=NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm")
    stop("\n'pls' must be an object of class 'plspm'")
  if (class(treemox) != "treemox")
    stop("\n'treemox' must be an object of class 'treemox'")
  if (nrow(pls$scores) != treemox$MOX$Size[1]) 
    stop("\n'pls' and 'treemox' have different number of observations")
  if (!is.null(X)) # if X available
  {
    if (is.null(pls$data))
    {
      if (!is_not_tabular(X))
        stop("\n'X' must be a numeric matrix or data frame")
      if (nrow(X) != nrow(pls$scores))
        stop("\n'pls' and 'X' have different number of rows")
    }
  } else { # if no X
    if (is.null(pls$data)) 
      stop("\nSorry, 'X' is missing. No dataset available")
  }
  
  # =======================================================
  # inputs setting
  # =======================================================  
  IDM <- pls$model$IDM
  blocks <- pls$model$blocks
  modes <- pls$model$specs$modes
  scheme <- pls$model$specs$scheme
  scaled <- pls$model$specs$scaled
  tol <- pls$model$specs$tol
  iter <- pls$model$specs$iter
  # data matrix DT
  if (!is.null(pls$data)) {
    DT <- pls$data
  } else {         
    # building data matrix 'DT'
    DT = get_manifests(X, blocks)
  }
  MOX <- treemox$MOX
  list.nodes <- treemox$list.nodes
  # identifying terminal nodes
  term.nodes <- which(MOX$Terminal=="yes") - 1    
  # number of terminal nodes
  tn <- length(term.nodes)
  # labels for terminal nodes
  tn.labs <- paste("Node", MOX$Node[term.nodes+1], sep="_")
  lvs <- ncol(IDM)
  obs <- nrow(DT)
  mvs <- pls$model$gens$mvs
  mvs.names <- pls$model$gens$mvs_names
  lvs.names <- pls$model$gens$lvs_names
  endo <- rowSums(IDM)
  endo[endo!=0] <- 1
  
  # =======================================================
  # PLS-PM of terminal nodes
  # =======================================================
  # prepare ingredients to store results
  TOW <- matrix(NA, mvs, 1+tn)  # TOW: table of outer weights
  TL <- matrix(NA, mvs, 1+tn)   # TL: table of loadings
  TR2 <- matrix(NA, lvs, 1+tn)  # TR2: table of R2
  # obtaining labels for the relations among latent variables
  path.labs <- NULL
  for (j in 1:lvs) {
    for (i in j:lvs) {
      if (IDM[i,j]==1) 
        path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
    }
  }
  # TPC: table of path coefficients
  TPC = matrix(NA, length(path.labs), 1+tn)
  
  ## model for root node
  TOW[,1] <- pls$outer_model$weight
  TL[,1] <- round(pls$outer_model$loading, 3)
  betas <- pls$path_coefs
  lvs.paths <- NULL
  for (j in 1:lvs)
    for (i in j:lvs)
      if (IDM[i,j]==1)
        lvs.paths <- c(lvs.paths, round(betas[i,j],3))
  TPC[,1] <- lvs.paths
  TR2[,1] <- round(pls$inner_summary$R2, 3)
  
  # calculating plspm results for terminal nodes
  for (k in 1:tn)
  {
    # elements of k-th terminal node
    elems.tn <- list.nodes[[term.nodes[k]]]
    DTN <- DT[elems.tn,]
    pls.tn <- get_pls_basic(DTN, IDM, blocks, pls$model$specs)
    # adding outer weights
    TOW[,k+1] <- pls.tn$out.weights
    # adding loadings
    TL[,k+1] <- pls.tn$loadings
    # obtaining relations among latent variables
    lvs.efs <- NULL
    betas <- pls.tn$path.coefs
    for (j in 1:lvs) {
      for (i in j:lvs) {
        if (IDM[i,j] == 1)
          lvs.efs <- c(lvs.efs, betas[i,j])
      }
    }
    TPC[,k+1] <- lvs.efs
    TR2[,k+1] <- pls.tn$R2
  }
  
  # =======================================================
  # Results
  # =======================================================
  # add names
  dimnames(TOW) <- list(mvs.names, c("Root_Node",tn.labs))
  dimnames(TL) <- list(mvs.names, c("Root_Node",tn.labs))
  dimnames(TPC) <- list(path.labs, c("Root_Node",tn.labs))
  dimnames(TR2) <- list(lvs.names, c("Root_Node",tn.labs))
  resul <- list(weights = TOW, 
                loadings = TL, 
                paths = TPC, 
                r2 = TR2, 
                IDM = IDM)
  class(resul) <- "treemox.pls"
  resul
}
