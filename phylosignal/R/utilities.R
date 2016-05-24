

#' Traits Values Simulation
#'
#' This function simulates the evolution of continuous characters along a phylogeny.
#' Traits values can be more or less influenced by the provided model.
#'
#' @param tree an object of class "\code{phylo}.
#' @param model the model to use ("\code{BM}" or "\code{OU}").
#' @param weight a numeric vector with values ranging between 0 and 1
#' giving the balance between pure random evolution and the selected model.
#' @param as.p4d logical. Should the phylogenetic tree and the simulated data coerced to a \code{phylo4d} object?
#' 
#' @return An object of class phyloSimSignal.
rTraitContWeight <- function(tree, model="BM", weight=1, as.p4d=FALSE){
  weight <- t(weight)
  trait.bm <- rTraitCont(tree, model=model)
  trait.rand <- sample(as.vector(trait.bm))
  names(trait.rand) <- names(trait.bm)
  trait.bm <- as.matrix(trait.bm)
  trait.rand <- as.matrix(trait.rand)
  trait.weight <- (trait.bm %*% weight) + (trait.rand %*% (1 - weight))
  colnames(trait.weight) <- paste("W", as.character(weight), sep="_")
  if(as.p4d){
    return(phylo4d(tree, tip.data=trait.weight, check.node.labels = "drop"))
  } else {
    return(trait.weight)
  }
}


#' phyloSimSignal subsetting
#'
#' This function subset all the components of a phyloSimSignal object following the methods and traits provided.
#'
#' @param x an object of class \code{phyloSimSignal} to subset.
#' @param methods a character vector giving the methods to keep.
#' @param traits a character vector giving the traits to keep.
#' 
#' @return A \code{phyloSimSignal} object
subsetPhyloSimSignal <- function(x, methods, traits){
  if(is.null(methods)){
    methods <- x$phylosim$stat.names
  }
  if(is.null(traits)){
    traits <- rownames(x$signal.mean)
  }
  methods.idx <- match(methods, x$phylosim$stat.names)
  traits.idx <- match(traits, rownames(x$signal.mean))
  x$signal.mean <- x$signal.mean[traits.idx, methods.idx]
  x$signal.qmin <- x$signal.qmin[traits.idx, methods.idx]
  x$signal.qmax <- x$signal.qmax[traits.idx, methods.idx]
  x$stat.signal$stat <- x$stat.signal$stat[traits.idx, methods.idx]
  x$stat.signal$pvalue <- x$stat.signal$pvalue[traits.idx, methods.idx]
  x$phylosim$sim.stat <- x$phylosim$sim.stat[ , methods.idx, ]
  x$phylosim$sim.pvalue <- x$phylosim$sim.pvalue[ , methods.idx, ]
  x$phylosim$stat.names <- x$phylosim$stat.names[methods.idx]
  
  return(x)
}


#' ID of direct descents
#'
#' This function return the node number (ID) of the direct descents of a given node
#'
#' @param phy an object of class "\code{phylo}.
#' @param node the ID of the node.
#' 
#' @return A vector of node numbers.
directDescents <- function(phy, node){
  res <- phy$edge[which(phy$edge[, 1] == node), 2]
  return(res)
}

#' Names of descents
#'
#' This function return all the names of the descents (tips) of a given node.
#'
#' @param phy an object of class "\code{phylo}.
#' @param node the ID of the node.
#' 
#' @return A vector tip labels.
descentsNames <- function(phy, node){
  tips <- phy$tip.label
  n.tips <- length(tips)
  if(node <= n.tips){
    res <- tips[node]
  } else {
    res <- extract.clade(phy, node = node)$tip.label
  }
  return(res)
}



#' Pairwise Distance from Regularly Distributed Points
#'
#' This function return a distance (or proximity) matrix between n points regularly distributed in 1 dimension.
#'
#' @param n the number of points
#' @param prox a logical indicating whether to return a matrix of proximity.
#' Default to \code{FALSE} so the function returns the matrix of distance.
#' 
#' @return A squared matrix. A matrix of class \code{dist} if \code{prox} is set to \code{FALSE}.
#' @examples
#' x <- distEq(5)
distEq <- function(n, prox = FALSE){
  M <- sapply(1:n, function(x) abs((1:n) - x))
  if(prox){
    M <- 1 / M
    diag(M) <- 0
  } else {
    M <- as.dist(M)
  }
  return(M)
}

#' Match matrix row/col names with phylo4d tips/traits
#'
#' Check for constitency of rows and columns names between a matrix of traits values and phylo4d object.
#'
#' @param x a matrix of data.
#' @param p4d a phylo4d object.
#' @param p4d.tips tips labels (relevant if p4d is \code{NULL}). 
#' @param p4d.traits traits labels (relevant if p4d is \code{NULL}).
#' @param subset a logical. Should the data matrix be subsetted using tips and traits labels.
#' 
#' @return The data matrix (eventually subsetted).
#' An error if no consistency between the data and the tree.
matchTipsAndTraits <- function(x, p4d = NULL, p4d.tips = NULL, p4d.traits = NULL, subset = TRUE){

  if(!is.null(p4d) & is(p4d, "phylo4d")){
    p4d.tips <- tipLabels(p4d)
    p4d.traits <- colnames(tdata(p4d))
  }
  
  if(!all(p4d.tips %in% rownames(x))){
    stop("Rows names of !!! do not match with tips names")
  }
  
  if(!all(p4d.traits %in% colnames(x))){
    stop("Columns names of !!! do not match with traits names")
  } 
  
  if(subset){
    x <- x[p4d.tips, p4d.traits]
  }
  return(x)
}


#' Set layout for plots
#' 
#' @param n.traits the number of traits in the layout
#' @param show.tip a logical indicating whether tip names are included in the layout
#' 
#' @rdname layouterize
.layouterize <- function(n.traits, show.tip){
  if(show.tip){
    res <- matrix(c(n.traits + 2, 1:(n.traits + 1)), nrow = 1)
  } else {
    res <- matrix(c(n.traits + 1, 1:(n.traits)), nrow = 1)
  }
  return(res)
}

#' Set layout for plots
#' 
#' @param tree.ratio the ratio of phylogenetic tree included in the layout
#' @param n.traits the number of traits in the layout
#' @param show.tip a logical indicating whether tip names are included in the layout
#' 
#' @rdname layouterizeRatio
.layouterizeRatio <- function(tree.ratio, n.traits, show.tip){
  if(!is.null(tree.ratio)){
    if(show.tip){
      res <- c(tree.ratio, rep((1 - tree.ratio) / (n.traits + 1), n.traits + 1))
    } else {
      res <- c(tree.ratio, rep((1 - tree.ratio) / n.traits, n.traits))
    }
  } else {
    if(show.tip){
      res <- rep(1, n.traits + 1)
    } else {
      res <- rep(1, n.traits)
    }      
  }
  return(res)
}

#' Color Palette
#' 
#' A simple color palette for gridplots.
#' 
#' @param n the number of colors to be in the palette.
white2red <- colorRampPalette(c("white", "red"))


#' Reordering vector or matrix of settings.
#'
#' Internal function of multiplot.phylo4d.
#' Reordering vector or matrix of settings and check for names consistencies.
#'
#' @param x a vector or a matrix to order.
#' @param n.tips number of tips.
#' @param n.traits number of traits
#' @param new.order a numeric vector giving the new order.
#' @param tips a character vector giving the tips labels (with the new order)
#' @param default the default value
#' 
#' @return An ordered vector or matrix.
#' An error if problem of consistency.
#' @rdname orderGrArg
.orderGrArg <- function(x, n.tips, n.traits, new.order, tips, default){
  x.dep <- deparse(substitute(x))
  if(is.vector(x)){
    if(is.null(names(x))){
      x <- rep(x, length.out = n.tips)
      x <- x[new.order]
    } else {
      if(any(tips %in% names(x))){
        y <- rep(default, n.tips)
        names(y) <- tips
        y[names(x)] <- x[names(x)]
        x <- y[tips]
      } else {
        stop(paste("Phylogenetic tips do not match with names of", x.dep))
      }
    }
  }
  if(is.matrix(x)){
    if(is.null(rownames(x))){
      x <- x[new.order, ]
    } else {
      if(any(tips %in% rownames(x))){
        y <- matrix(default, nrow = nrow(x), ncol = ncol(x))
        rownames(y) <- tips
        y[rownames(x),] <- x[rownames(x),]
        x <- y[tips, ]
      } else {
        stop(paste("Phylogenetic tips do not match with row names of", x.dep))
      }
    }
  }
  x <- matrix(x, nrow = n.tips, ncol = n.traits)
  return(x)
}



# .simpleCat <- function(x){
#   if(!is.null(names(x))){
#     x.names <- names(x)
#   } else {
#     x.names <- NULL
#   }
#   x <- as.factor(x)
#   x <- as.numeric(x)
#   names(x) <- x.names
#   return(x)
# }

#' Palette of evenly distributed colors
#' 
#' This function generates a vector of n colors evenly distributed in the RGB space.
#' Usefull to create palettes of distinct colors.
#' 
#' @param n the number of colors to be in the palette.
#' 
#' @return a vector of hexadecimal colors.
#' @export
evenColors <- function(n){
  rgb.cube <- expand.grid(seq(0, 255, by = 15),
                          seq(0, 255, by = 15),
                          seq(0, 255, by = 15))
  rgb.km <- kmeans(rgb.cube, n)
  res <- mapply(rgb, rgb.km$centers[, 1], rgb.km$centers[, 2], rgb.km$centers[, 3],
                MoreArgs = list(maxColorValue = 255))
  return(res)
}
