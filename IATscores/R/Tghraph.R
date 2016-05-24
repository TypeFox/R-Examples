Tgraph <- function(mcmp, alpha = .05, horizorder = NULL)
{
  x <- mcmp$Analysis
  # from and to
  nd <- row.names(x)
  nd <- sapply(nd, str_split, pattern = "-")
  from <- unname(sapply(nd, str_trim)[1,])
  to <- unname(sapply(nd, str_trim)[2,])
  edgelist <- data.frame(cbind(from, to, x))
  edgelist$from <- as.character(edgelist$from)
  edgelist$to <- as.character(edgelist$to)
  
  # keep only significant edeges
  edgelist <- edgelist[edgelist$p.Value < alpha,]
  
  # if the estimated value is negative, just reverse the direction of the arrow
  tmp <- edgelist[edgelist[, "Estimator"] < 0, "to"]
  edgelist[edgelist[, "Estimator"] < 0, "to"] <-
    edgelist[edgelist[, "Estimator"] < 0, "from"]
  edgelist[edgelist[, "Estimator"] < 0, "from"] <- tmp
  edgelist <- edgelist[, c("from", "to", "Estimator")]
  edgelist$Estimator <- abs(edgelist$Estimator)
  
  # from edgelist to adjacency / weight matix
  wmat <- getWmat(edgelist)
  amat <- ceiling(wmat)
  amat <- transitive.reduction(amat)
  
  # FIND HIERARCHICAL LEVELS IN THE GRAPH
  lev <- hlevels(amat)
  
  # build the hierarchical layout for qgraph
  if(length(lev) > 1)
    layout <- lev2layout(lev)
  wmat[amat == 0] <- 0
  
  if(!is.null(horizorder))
  {
    # if a preferred horizontal order is specified, it is applied
    # 1. establish a vertical order
    layout <- data.frame(layout)
    names(layout) <- c("x", "y")
    layout[, "names"] <- rownames(amat)
    # 2. within levels of the vertical order, reorder the rows
    for(i in unique(layout$y))
    {
      xs <- unique(layout[layout$y == i, "x"])
      nds <- unique(layout[layout$y == i, "names"])
      nds2 <- match(nds, horizorder[horizorder %in% nds])
      nds2 <- nds2[!is.na(nds2)]
      newxs <- xs[nds2]
      layout[layout$y == i, "x"] <- newxs
    }
    layout <- as.matrix(layout[,c("x", "y")])
  }
  
  list("wmat" = wmat, "amat" = amat, "layout" = layout)
}


hlevels <- function(amat)
{
  # starting from an adjacency matrix representing a DAG, this function
  # decides the layers to plot the nodes. The first layer are the nodes
  # that do not receive edges. The second layer are those nodes that rececive
  # from the first layer. The thir layer recceives from the second and so on.
  
  # amat = directed matrix with hierarchical layout (e.g., DAG)
  # after applying  nem::transitive.reduction(amat)
  
  # remove nodes that do not send or receive arrows
  toremove <- colSums(amat) == 0 & rowSums(amat) == 0
  amat <- amat[!toremove, !toremove]
  
  # lev is a list that includes the nodes for each level
  lev <- list()
  nodes <- (1:ncol(amat))
  
  # those that are not dominated are in the upper (1st) level
  lev[[1]] <- nodes[colSums(amat) == 0]
  cont <- 2
  
  while(length(unlist(lev)) != ncol(amat))
  {
    if(length(lev[[cont-1]]) == 1)
    {
      index <- amat[lev[[cont-1]],]
    }  else {
      index <- colSums(amat[lev[[cont-1]],] != 0)
    }
    
    if(sum(index != 0) > 1)
    {
      # if any of the nodes in a level receives incoming arrows from the others
      # of the same level, it goes in the subsequent level
      nextlev <- colSums(amat[as.logical(index), as.logical(index)])
      if(any(nextlev != 0)) index[index!=0][as.logical(nextlev) != 0] <- 0
    }
    
    lev[[cont]] <- nodes[as.logical(index)]
    cont <- cont+1
  }
  if(length(lev) <= 1) stop("no significant contrast in input")
  lev
}


lev2layout <- function(lev)
{
  # creates a qgraph layout from a series of levels.
  layout <- matrix(0, ncol = 2, nrow = length(unlist(lev)))
  stepy <- 2/(length(lev)-1)
  ys <- seq(from = 1, to = -1, by = -stepy)
  
  for(i in 1:length(lev))
  {
    # y axis represents the level
    layout[lev[[i]], 2] <- ys[i]
    
    # x axis = equally spaced
    if(length(lev[[i]]) == 1)
    {
      layout[lev[[i]], 1] <- 0
    } else {
      stepx <- 2/(length(lev[[i]]) - 1)
      xs <- seq(from = -1, to = 1, by = stepx)
      layout[lev[[i]], 1] <- xs
    }
  }
  layout
}