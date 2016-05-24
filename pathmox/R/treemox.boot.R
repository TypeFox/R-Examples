#' @title Bootstrapping validation for PATHMOX or TECHMOX trees
#' 
#' @description
#' Performs bootstrapping validation on path coefficients of terminal nodes from
#' a PATHMOX or TECHMOX tree
#' 
#' @details
#' The default number of re-samples is 100. However, \code{br} can be specified
#' in a range from 50 to 500. \cr When the object \code{pls} does not contain a
#' data matrix (i.e. \code{pls$data=NULL}), the user must provide the data
#' matrix or data frame in \code{X}.
#' 
#' @param pls An object of class \code{"plspm"} returned by \code{\link{plspm}}.
#' @param treemox An object of class \code{"treemox"} returned by either
#' \code{\link{pathmox}} or \code{\link{techmox}}.
#' @param X Optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{pls}.
#' @param br An integer indicating the number bootstrap resamples (\code{br=100}
#' by default).
#' @return An object of class \code{"bootnodes"}. Basically a list with the
#' following results:
#' @return \item{PC}{Matrix of original path coefficients for the root node and
#' the terminal nodes.}
#' @return \item{PMB}{Matrix of bootstrap path coefficients (mean value) 
#' for the root node and the terminal nodes.}
#' @return \item{PSB}{Matrix of bootstrap standard errors of path coefficients 
#' for the root node and the terminal nodes.}
#' @return \item{PP05}{Matrix of 0.05 bootstrap percentile of path coefficients 
#' for the root node and the terminal nodes.}
#' @return \item{PP95}{Matrix of 0.95 bootstrap percentile of path coefficients 
#' for the root node and the terminal nodes.}
#' @author Gaston Sanchez
#' @seealso \code{\link{pathmox}}, \code{\link{techmox}},
#' \code{\link{treemox.pls}}.
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
#'  # re-ordering those segmentation variables with ordinal scale 
#'  # (Age and Education)
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
#'  mob_nodes_boot = treemox.boot(mob_pls, mob_pathmox)
#'
#'  # plot of results for path coefficient number 12
#'  plot(mob_nodes_boot, pc=12)  
#'  }
#'
treemox.boot <- function(pls, treemox, X=NULL, br=100)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm")
    stop("\n'pls' must be an object of class 'plspm'")
  if (class(treemox) != "treemox")
    stop("\n'treemox' must be an object of class 'treemox'")
  if (nrow(pls$scores) != treemox$MOX$Size[1]) 
    stop("\n'pls' and 'pathmox' have different number of observations")
  if (!is.null(X)) # if X available
  {
    if (is.null(pls$data))
    {
      if (!is.matrix(X) && !is.data.frame(X))
        stop("\n'X' must be a numeric matrix or data frame.")
      if (nrow(X) != nrow(pls$scores))
        stop("\n'pls' and 'X' have different number of rows.")
      if (nrow(X) != pathmox$MOX$Size[1])
        stop("\n'X' and 'pathmox' have different number of observations")
    }
  } else { # if no X
    if (is.null(pls$data))
      stop("\nSorry, 'X' is missing. No dataset available.")
  }
  if (mode(br)!="numeric" || length(br)!=1 || (br%%1)!=0 ||
    br<50 || br>500) {
    warning("\ninvalid 'br'; default 'br=100' is used.")   
    br <- 100
  } 
  
  # =======================================================
  # inputs setting
  # =======================================================  
  IDM <- pls$model$IDM
  blocks <- pls$model$blocks
  scaled <- pls$model$specs$scaled

  # data matrix DT
  if (!is.null(pls$data)) {
    DT <- pls$data
  } else {
    # building data matrix 'DT'
    DT = get_manifests(X, blocks)
  }
  lvs <- ncol(IDM)
  obs <- nrow(DT)
  endo <- rowSums(IDM)
  endo[endo!=0] <- 1
  obs.names = rownames(DT)
  mvs.names = colnames(DT)
  lvs.names = colnames(IDM)
  
  # =======================================================
  # Parameters from pathmox
  # =======================================================  
  ### Parameters from pathmox
  MOX <- treemox$MOX
  list.nodes <- treemox$list.nodes
  term.nodes = which(MOX$Terminal=="yes") - 1# identifying terminal nodes    
  tn = length(term.nodes)# number of terminal nodes
  tn.labs = paste("Node", MOX$Node[term.nodes+1], sep="_")# labels for terminal 
  
  path.labs <- NULL
  for (j in 1:lvs) {
    for (i in j:lvs) {
      if (IDM[i,j]==1) 
        path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep="")) 
    }
  }
  TPC <- matrix(NA, length(path.labs), 1+tn)# TPC: table of path coefficients   
  
  ### Path coefficients for Root Node
  betas <- pls$path_coefs
  lvs.paths <- NULL
  for (j in 1:lvs) {
    for (i in j:lvs) {
      if (IDM[i,j]==1)
        lvs.paths <- c(lvs.paths, round(betas[i,j],3))
    }
  }
  TPC[,1] <- lvs.paths
  
  # calculating plspm results for terminal nodes
  for (k in 1:tn)
  {        
    elems.tn <- list.nodes[[term.nodes[k]]]# elements of k-th terminal node
    DTN <- DT[elems.tn,]
    pls.tn <- get_pls_basic(DTN, IDM, blocks, pls$model$specs)
    # obtaining relations among latent variables
    lvs.efs <- NULL
    betas <- pls.tn$path.coefs
    for (j in 1:lvs) {
      for (i in j:lvs) {
        if (IDM[i,j]==1)
          lvs.efs <- c(lvs.efs, round(betas[i,j],3))
      }
    }
    TPC[,k+1] <- lvs.efs
  }
  dimnames(TPC) <- list(path.labs, c("Root_Node",tn.labs))
  
  # =======================================================
  # Bootstrapping
  # =======================================================  
  # default number of samples br=100
  bootnum <- br
  Path.meanboot = matrix(NA, sum(IDM), tn+1)
  Path.steboot = matrix(NA, sum(IDM), tn+1)
  Path.perc05 = matrix(NA, sum(IDM), tn+1)
  Path.perc95 = matrix(NA, sum(IDM), tn+1)
  
  # =============== Bootstrapping for Root Node ===========
  PATHS <- matrix(NA, bootnum, sum(IDM))
  i <- 1
  while (i <= bootnum)
  {
    boot.obs <- sample.int(nrow(DT), size=nrow(DT), replace=TRUE)
    DM.boot <- DT[boot.obs,]
    # scaling boot sample
    if (scaled) {
      sd.XB <- sqrt((nrow(DM.boot)-1)/nrow(DM.boot)) * apply(DM.boot, 2, sd)
      X.boot <- scale(DM.boot, scale=sd.XB)
    } else {
      X.boot <- scale(DM.boot, scale=FALSE)
    }
    colnames(X.boot) <- mvs.names
    # calculating boot model parameters 
    w.boot <- get_weights(X.boot, IDM, blocks, pls$model$specs)
    if (is.null(w.boot)) {
      i <- i - 1
      next
    }
    Y.boot <- X.boot %*% w.boot[[2]]
    pathmod <- get_paths(IDM, Y.boot)
    P.boot <- pathmod[[2]]
    PATHS[i,] <- as.vector(P.boot[which(IDM==1)])
    i <- i + 1
  }
  Path.meanboot[,1] <- apply(PATHS, 2, mean)
  Path.steboot[,1] <- apply(PATHS, 2, sd)
  Path.perc05[,1] <- apply(PATHS, 2, function(x) quantile(x,p=.05))
  Path.perc95[,1] <- apply(PATHS, 2, function(x) quantile(x,p=.95))
  
  # =======================================================
  # Bootstrapping of terminal nodes
  # =======================================================  
  for (k in 1:tn)
  {        
    elems.tn = list.nodes[[term.nodes[k]]]# elements of k-th terminal node
    DT.new = DT[elems.tn,]
    PATHS = matrix(NA, bootnum, sum(IDM))
    i <- 1
    while (i <= bootnum)
    {            
      boot.obs <- sample.int(nrow(DT.new), size=nrow(DT.new), replace=TRUE)
      DM.boot <- DT.new[boot.obs,]
      # scaling boot sample
      if (scaled) {
        sd.XB <- sqrt((nrow(DM.boot)-1)/nrow(DM.boot)) * apply(DM.boot, 2, sd)
        X.boot <- scale(DM.boot, scale=sd.XB)
      } else {
        X.boot <- scale(DM.boot, scale=FALSE)
      }
      colnames(X.boot) <- mvs.names
      # calculating boot model parameters 
      w.boot <- get_weights(X.boot, IDM, blocks, pls$model$specs)
      if (is.null(w.boot)) {
        i <- i - 1
        next
      }
      Y.boot <- X.boot %*% w.boot[[2]]
      pathmod <- get_paths(IDM, Y.boot)
      P.boot <- pathmod[[2]]
      PATHS[i,] <- as.vector(P.boot[which(IDM==1)])
      i <- i + 1
    }
    Path.meanboot[,(k+1)] = apply(PATHS, 2, mean)
    Path.steboot[,(k+1)] = apply(PATHS, 2, sd)
    Path.perc05[,(k+1)] = apply(PATHS, 2, function(x) quantile(x,p=.05))
    Path.perc95[,(k+1)] = apply(PATHS, 2, function(x) quantile(x,p=.95))
  } # end for 'tn' terminal nodes
  
  # =======================================================
  # Results
  # =======================================================
  # add names
  dimnames(Path.meanboot) <- list(path.labs, c("Root_Node",tn.labs))
  dimnames(Path.steboot) <- list(path.labs, c("Root_Node",tn.labs))
  dimnames(Path.perc05) <- list(path.labs, c("Root_Node",tn.labs))
  dimnames(Path.perc95) <- list(path.labs, c("Root_Node",tn.labs))
  # results
  bootnodes <- list(PC = TPC, 
                    PMB = Path.meanboot, 
                    PSB = Path.steboot, 
                    PP05 = Path.perc05, 
                    PP95 = Path.perc95)
  class(bootnodes) <- "bootnodes"
  bootnodes
}

