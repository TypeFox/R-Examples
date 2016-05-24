#' @title TECHMOX Algorithm: Segmentation Trees in PLS Path Modeling
#' 
#' @description
#' The function \code{techmox} calculates a binary segmentation tree following
#' the TECHMOX algorithm.  In contrast, \code{fix.techmox} obtains a supervised
#' TECHMOX tree in the sense of allowing the user to interactively fix the
#' partitions along the construction process of the tree. \cr
#' 
#' @details
#' The argument \code{EXEV} must be a data frame containing segmentation
#' variables as factors (see \code{\link{factor}}). The number of rows in
#' \code{EXEV} must be the same as the number of rows in the data used in
#' \code{pls}. \cr
#' 
#' The argument \code{size} can be defined as a decimal value (i.e. proportion
#' of elements inside a node), or as an integer (i.e. number of elements inside
#' a node). \cr
#' 
#' When the object \code{pls} does not contain a data matrix (i.e.
#' \code{pls$data=NULL}), the user must provide the data matrix or data frame 
#' in \code{X}.
#' 
#' @aliases techmox fix.techmox
#' @param pls An object of class \code{"plspm"} returned by \code{\link{plspm}}.
#' @param EXEV A data frame of factors contaning the segmentation variables.
#' @param X Optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{pls}.
#' @param signif A numeric value indicating the significance threshold of the
#' F-statistic. Must be a decimal number between 0 and 1.
#' @param size A numeric value indicating the minimum size of elements inside a
#' node.
#' @param deep An integer indicating the depth level of the tree. Must be an
#' integer greater than 1.
#' @param tree A logical value indicating if the tree should be displayed
#' (\code{TRUE} by default).
#' @return An object of class \code{"treemox"}. Basically a list with the
#' following results:
#' @return \item{MOX}{Data frame containing the results of the segmentation tree.}
#' @return \item{FT}{Data frame containing the results of the F-tests for each node
#' partition.}
#' @return \item{candidates}{List of data frames containing the candidate splits of
#' each node partition.}
#' @return \item{list.nodes}{List of elements for each node.}
#' 
#' @author Gaston Sanchez
#' @references Sanchez, G. (2009) \emph{PATHMOX Approach: Segmentation Trees in
#' Partial Least Squares Path Modeling.} Doctoral Dissertation. 
#' 
#' \url{http://www.gastonsanchez.com/Pathmox_Approach_Thesis_Gaston_Sanchez.pdf}
#' @seealso \code{\link{pathmox}}, \code{\link{plot.treemox}}
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
#'  # Techmox Analysis
#'  mob_techmox = techmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)
#'  }
#'
techmox <-
function(pls, EXEV, X=NULL, signif=.05, size=.10, deep=2, tree=TRUE)
{ 
  # =======================================================
  # checking arguments
  # =======================================================
  valid = check_pathmox_args(pls = pls, EXEV = EXEV, X = X, 
                             signif = signif, size = size, deep = deep, tree = tree)
  signif = valid$signif
  size = valid$size
  deep = valid$deep
  tree = valid$tree
  
  # =======================================================
  # inputs setting
  # =======================================================  
  N = nrow(EXEV)
  IDM = pls$model$IDM
  blocks = pls$model$blocks

  # data matrix DT
  if (!is.null(pls$data)) {
    DT <- pls$data
  } else {
    # building data matrix 'DT'
    DT = get_manifests(X, blocks)
  }
  
  list.nodes <- NULL
  FT.global <- NULL
  FT.partial <- NULL
  FT.labs <- NULL
  FT.otras <- NULL
  type.exev <- rep(0,ncol(EXEV))
  treat.exev <- rep("binary", ncol(EXEV))
  for (i in 1:length(type.exev))
  {
    type.exev[i] <- ifelse(is.ordered(EXEV[,i]), "ord", "nom")
    if (nlevels(EXEV[,i]) > 2) {
      if (is.ordered(EXEV[,i])) {
        treat.exev[i] = "ordinal" 
      } else treat.exev[i]="nominal"
    }
  }
  df.exev <- data.frame(Nlevels = unlist(lapply(EXEV,nlevels)), 
                        Ordered = unlist(lapply(EXEV,is.ordered)), 
                        Treatment = treat.exev)
  
  # =======================================================
  # XEXELOA NELHUATL
  # =======================================================  
  elements <- 1:nrow(EXEV)   # nelhuatl elems
  elemnod <- rep(0, N)
  nv <- 0
  spli <- get_xexeloa(pls, DT, EXEV, type.exev, elemnod, nv, size, "techmox")
  nodes <- 2:3
  parent <- c(1,1)
  nodes.level <- rep(1,2)
  FT.global <- rbind(FT.global, spli$inner.global) # global
  FT.partial <- c(FT.partial, list(spli$inner.partial$F.test)) # partial
  FT.labs <- c(FT.labs, spli$inner.partial$exev)
  FT.otras <- c(FT.otras, list(spli$otras))# otras
  exvar <- rep(spli$inner.global[1], 2)
  best <- spli$inner.partial$p.val    
  if (deep==1 || best>signif) {
    type.node=c("leaf","leaf")
    terminal.node <- c("yes", "yes")
  } else {
    type.node=c("node","node")  
    terminal.node <- c("no", "no")
  }
  lef.split <- spli$list.elems[[1]]# a.cuauh
  rig.split <- spli$list.elems[[2]]# b.cuauh
  list.nodes <- list(lef.split, rig.split)
  size.node <- c(length(lef.split), length(rig.split))
  elemnod[lef.split] <- 2# label elems a.cuauh
  elemnod[rig.split] <- 3# label elems b.cuauh   
  categories <- list(spli$categs[[1]], spli$categs[[2]])
  
  ### nelhuatl xexeloa
  if (size<1)  size.limit = ceiling(size*N)  else  size.limit = size
  nivel <- 1   
  ### cuahuitl ohtli
  repeat
  {
    if (deep==1 || best>signif)
      break
    num.nods.niv <- 2^nivel   # nodos q ai n teoria
    theo.nods.niv <- (2^nivel):(2^(nivel+1)-1)  # nodos q ai n teoria
    nods.niv <- intersect(nodes, theo.nods.niv)   # nodos q ai
    for (nd in 1:length(nods.niv))   # p/c cuauhmaitl
    {
      nv <- nods.niv[nd]   # su etiketado
      indivs <- which(elemnod == nv)
      # cuahuitl aturat 
      n.size <- length(indivs)
      if (n.size <= size.limit)# aturat
      {
        type.node[which(nodes == nv)] <- "leaf"    
        terminal.node[which(nodes == nv)] <- "yes"    
        next  
      } else # segueix
      {
        DT.set <- DT[indivs,]
        x.pls <- get_pls_basic(DT.set, IDM, blocks, pls$model$specs)
        # xexeloa
        spli <- get_xexeloa(x.pls, DT.set, EXEV, type.exev, elemnod, 
                            nv, size, "techmox")    
        best <- spli$inner.partial$p.val    
        # ai gallo o nel?
        if (best > signif)   # aturat
        {
          type.node[which(nodes==nv)] <- "leaf"    
          terminal.node[which(nodes==nv)] <- "yes"    
        } else   
        {
          FT.global <- rbind(FT.global, spli$inner.global)
          FT.partial <- c(FT.partial, list(spli$inner.partial$F.test))
          FT.labs <- c(FT.labs, spli$inner.partial$exev)
          FT.otras <- c(FT.otras, list(spli$otras))
          exvar <- c(exvar, rep(spli$inner.global[1],2))
          nodes <- c(nodes, 2*nv, 2*nv+1)   
          parent <- c(parent, rep(nv,2))
          nodes.level <- c(nodes.level, rep(nivel+1,2))
          type.node <- c(type.node, rep("node",2))   
          terminal.node <- c(terminal.node, rep("no",2))  
          lef.split <- spli$list.elems[[1]]
          rig.split <- spli$list.elems[[2]]
          ultim <- length(type.node)
          if (length(lef.split) <= size.limit)
          {
            type.node[ultim-1] <- "leaf"
            terminal.node[ultim-1] <- "yes"
          }
          if (length(rig.split) <= size.limit)
          {
            type.node[ultim] <- "leaf"
            terminal.node[ultim] <- "yes"
          }
          list.nodes <- c(list.nodes, list(lef.split, rig.split))
          siz.nd <- c(length(lef.split), length(rig.split))
          size.node <- c(size.node, siz.nd)
          elemnod[lef.split] <- 2*nv
          elemnod[rig.split] <- 2*nv+1
          cat <- length(categories)
          categories[[cat+1]] <- spli$categs[[1]] 
          categories[[cat+2]] <- spli$categs[[2]] 
        }
      }
    }
    nivel <- nivel + 1 
    if (nivel >= deep) 
      break
    nods.nxt.niv <- (2^nivel):(2^(nivel+1)-1)# nodos q ai n teoria
    nnn <- which(nodes %in% nods.nxt.niv)  # (nnn: nods.nxt.niv)
    if (length(nnn)==0)   # aturat cullons
      break
  }
  # end repeat
  
  # cuauhmaitl terminal ultim xexeloa 
  nods.nxt.niv <- (2^nivel):(2^(nivel+1)-1)
  nnn <- which(nodes %in% nods.nxt.niv)
  if (length(nnn) > 0)
    terminal.node[nnn] <- rep("yes", length(nnn))    
  excat <- NULL
  for (h1 in 1:length(categories))
    excat <- rbind(
      excat, 
      paste(as.vector(categories[[h1]]), sep="", collapse="/"))    
  
  # =======================================================
  # Results
  # =======================================================  
  names(list.nodes) <- nodes
  nodes <- c(1, nodes)
  parent <- c(0, parent)
  nodes.level <- c(0, nodes.level)
  type.node <- c("root", type.node)
  terminal.node <- c("no", terminal.node)
  size.node <- c(N, size.node)
  Perc <- 100*(size.node/N)    
  Variable <- c(NA, as.vector(unlist(exvar)))
  Categ <- c(NA, excat)
  model <- list(mox = "techmox", 
                signif = signif, 
                size = size, 
                deep = deep, 
                df.exev = df.exev)
  node.nums <- unique(parent)[-c(1,2)]
  node.labs <- c("Root", paste(rep("Node",length(node.nums)),node.nums,sep="") )   
  names(FT.otras) <- node.labs
  names(FT.partial) <- paste(node.labs, FT.labs, sep="_")
  MOX <- data.frame(Node = nodes, 
                    Parent = parent, 
                    Depth = nodes.level, 
                    Type = type.node, 
                    Terminal = terminal.node, 
                    Size = size.node, 
                    Percent = Perc, 
                    Variable = Variable, 
                    Category = Categ )
  mox.resul <- list(MOX = MOX, 
                    FT = FT.partial, 
                    candidates = FT.otras, 
                    list.nodes = list.nodes, 
                    model = model)
  class(mox.resul) <- "treemox"
  if (tree) plot.treemox(mox.resul)
  mox.resul
}
