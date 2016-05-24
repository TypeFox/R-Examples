#' Plot recombination frequencies
#' 
#' Graphical representation of recombination frequencies to support supervised estimation
#' of the numbers of clusters
#'  
#' @param rf Matrix of pairwise recombination frequencies.  
#' @param plottype Default is "dendrogram". Any other value will plot the recombination frequencies.
#' @param method Default is "single", which is used for the hierarchical clustering.
#' @param cex.axis Size of axis labels in image plot.
#' @param ... arguments are forwarded to \code{image}.
#' @return None.
#' @import graphics
#' @examples
#' data(simTetra)
#' simTetraGen <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetraGen, 4)
#' plotRf(rfMat)
#' @export
plotRf <- function(rf, plottype = "dendrogram", method = "single", cex.axis = 1, ...){
  tree<-hclust(as.dist(rf), method = method)
  if(plottype == "dendrogram"){
    plot(tree, xlab = "Linkage groups", ...)
  }else{
    image(rf[tree$order, tree$order], xaxt = 'n', yaxt = 'n', ...)
    axis(side = 1, at=seq(0, 1, length.out = length(tree$order)),
         labels = rownames(rf)[tree$order], las = 2, cex.axis = cex.axis)
    axis(side = 2, at = seq(0, 1, length.out = length(tree$order)),
         labels = rownames(rf)[tree$order], las = 2, cex.axis = cex.axis)    
  }
}

#' Plotting one or two linkage maps
#' 
#' Visualization of one or two linkage maps.
#' Used as comparison between two different maps (e.g. different parameters or linkage mapping tools).
#'   
#' @param map1 Numeric vector with marker positions.
#' @param map2 Optional second map for comparison.
#' @param cex Font size in the figure.
#' @param labels Labels for the two blocks
#' @param ... arguments are forwarded to \code{plot}.
#' @return None. Plotting only.
#' @import graphics
#' @examples
#' data(simTetra)
#' simTetraGen <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetraGen, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' split <- sortLeafs(rfMat, split)
#' map <- pullMap(rfMat, split = split) 
#' plotChr(map[[1]])  
#' @export

plotChr <- function(map1, map2 = NULL, cex = 1, labels = c("Map 1", "Map 2"), ...){
  nMark1 <- length(map1)
  x1 = 0.9
  x2 = 1.1
  x3 = 1.9
  x4 = 2.1
  t1 = 0.5
  t2 = 2.5
  #plot single chromosome
  if(is.null(map2)){
    plot(0, type = "n", axes = FALSE, ann = FALSE, xlim = c(0, 3), 
         ylim = c(0, max(map1)), ylab = "cM", xlab = labels[1], ...)
    lines(c(x1, x1), c(min(map1), max(map1)), type = "l")
    lines(c(x2, x2), c(min(map1), max(map1)), type = "l")
    segments(x0 = rep(x1, nMark1), x1 = rep(x2, nMark1), y0 = map1, y1 = map1)
    text(x = t1, y = map1, labels = names(map1), cex = cex)
  }else{ #compare two maps
    plot(0, type = "n", axes = FALSE, xlim = c(0, 3), ylim = c(0, max(c(map1, map2))), 
         ylab = "cM", xlab = labels[1], ...)   
    lines(c(x1, x1), c(min(map1), max(map1)), type = "l")
    lines(c(x2, x2), c(min(map1), max(map1)), type = "l")
    segments(x0 = rep(x1, nMark1), x1 = rep(x2, nMark1), y0 = map1, y1 = map1)
    text(x = t1, y = map1, labels = names(map1), cex = cex)
    #second map
    lines(c(x3, x3), c(min(map2), max(map2)), type = "l")
    lines(c(x4, x4), c(min(map2), max(map2)), type = "l")
    nMark2 <- length(map2)
    segments(x0 = rep(x3, nMark2), x1 = rep(x4, nMark2), y0 = map2, y1 = map2) 
    text(x = t2, y = map2, labels = names(map2), cex = cex)
    map1InMap2 <- names(map1) %in% names(map2)
    nMarkBoth <- sum(map1InMap2)
    segments(x0 = rep(x2, nMarkBoth), x1 = rep(x3, nMarkBoth), y0 = map1[map1InMap2], 
             y1 = map2[names(map1)[map1InMap2]], col = 1 )
    
    axis(1, line = NA, at = 1:2, labels = labels, lwd = 0)
    axis(2, at = seq(0, max(c(map1, map2)), 10))
  }
}


#' Create a gray scale tanglegram
#' 
#' Create tanglegram. Removes markers, that are not in both trees.
#' Calculates alternating light and dark shades of grey.
#'     
#' @param dend1 First dendrogram. Required.
#' @param dend2 Second dendrogram. Required.
#' @param cutheight The height, at which dend1 is cut. Influences number of colors.
#' @param k Number of desired linkage groups.
#' @param ncol Number of desired colors.
#' @param ... Other parameters are forwarded to the tanglegram command.
#' @return None. Plotting only.
#' @examples
#' data(simTetra)
#' simTetraGen <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetraGen, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' split <- sortLeafs(rfMat, split)
#' map <- pullMap(rfMat, split = split)  
#' dend <- map2dend(map)  
#' maketangle(dend, dend, cutheight = 500, k = 7, ncol = 7)
#' @export
maketangle<-function(dend1, dend2, cutheight, k = NULL, ncol = k, ...){
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("dendextend needed for this function to work. Please install it.",
         call. = FALSE)
  }
  dendlist <- dendextend::intersect_trees(dend1, dend2)
  split <- dendextend::cutree(dendlist[[1]], h = cutheight)
  if(missing(k)) k <- max(split)
  dendextend::tanglegram(dendlist, color_lines = grDevices::gray.colors(ncol)[makealtord(k)][split], ...)
}


#' Creates vectors with highly distant neighbors
#' 
#' Creates a vector of numbers 1 to n, where the neighbors are as distant as possible.
#' For the grey scale plot, that guarantees, that the shades of grey are easy to distinguish.
#' For instance, 1, 4, 2, 5, 3; all numbers have a distance of at least 2 and where possible 3.     
#' @param n Length of vector.
#' @return Vector of length n.
#' @keywords internal
makealtord <- function(n = 3){  
  out<-rep(0, n)
  out[seq(1, n, 2)] <- 1:ceiling(n / 2)
  out[seq(2, n, 2)] <- (ceiling(n / 2) + 1):n
  return(out)
}
