
#' Phylogenetic correlogram
#' 
#' This function computes a phylogenetic correlogram.
#' 
#' @param p4d a \code{phylo4d} object.
#' @param trait the traits in the \code{phylo4d} object to use for the correlogram.
#' Can be a character vector giving the name of the traits or numbers giving the column index
#' in the table of the data slot of the \code{phylo4d} object.
#' @param dist.phylo a matrix of phylogenetic distances or a character string specifying a method to compute it.
#' See Details.
#' @param sigma a numeric value giving the standard deviation of the normal distribution used
#' to compute the matrix of phylogenetic weights. If \code{NULL} (default),
#' the function computes a value from the phylogeny.
#' @param n.points an integer giving the number of points at which to compute the correlogram's statistics.
#' @param ci.bs an integer giving the number of bootstrap replicates for confidence interval estimation.
#' @param ci.conf a value between 0 and 1 giving the confidence level of the confidence interval.
#' 
#' @details
#' This function computes a correlogram on a continuous scale of phylogenetic distance.
#' This is achieved by using a collection of specific phylogenetic weights matrices generated
#' with the "\code{lag-norm}" method of \code{\link{phyloWeights}} and different values of "\code{mu}".
#' 
#' The confidence envelope is computed by bootstrapping.
#' At each iteration, the autocorrelation is re-estimated after re-standardization of the
#' matrix of phylogenetic weights.
#' The nonparametric confidence intervals are computed at each lag by first order normal approximation.
#' Intervals are constrained between 0 and 1.
#' 
#' If there is one trait, the function computes Moran's I. If there is more than one trait,
#' the function computes the Mantel's statistic (Oden and Sokal 1986).
#' 
#' If "\code{dist.phylo}" is a character string,
#' the phylogenetic distance matrix is computed internally
#' using the function \code{\link[adephylo]{distTips}} from the package \pkg{adephylo}.
#' See \code{\link[adephylo]{distTips}} for details about the methods.
#' 
#' @return An object of class "\code{phylocorrelogram}".
#' 
#' @seealso \code{\link{plot.phylocorrelogram}},
#' \code{\link[ape]{correlogram.formula}} in \pkg{ape} for correlograms based on taxonomic levels.
#' 
#' @examples
#' \dontrun{
#' data(navic)
#' pc <- phyloCorrelogram(navic)
#' plot(pc)
#' }
#' @references 
#' Oden N.L. & Sokal R.R. (1986) Directional Autocorrelation: An Extension of Spatial Correlograms to Two Dimensions. Systematic Zoology 35, 608-617.
#' 
#' @export
phyloCorrelogram <- function(p4d, trait = names(tdata(p4d)),
                             dist.phylo = "patristic", sigma = NULL,
                             n.points = 100, ci.bs = 1000, ci.conf = 0.95){
  
  p4 <- extractTree(p4d)
  phy <- as(p4, "phylo")
  new.order <- phy$edge[, 2][!phy$edge[, 2] %in% phy$edge[, 1]]
  tips <- phy$tip.label[new.order]
  n.tips <- length(tips)
  X <- tdata(p4d, type = "tip")
  X <- X[tips, trait]
  X <- scale(X)
  X <- as.data.frame(X)
  colnames(X) <- trait
  n.traits <- ncol(X)
  if(is.numeric(trait)){
    trait <- names(tdata(p4d))[trait]
  }
  
  if(is.vector(dist.phylo) & is.character(dist.phylo)){
    dist.phylo <- match.arg(dist.phylo, c("patristic", "nNodes", "Abouheif", "sumDD"))
    D <- distTips(phy, method = dist.phylo)
    D <- as.matrix(D)
    D <- D[tips, tips]
  } else {
    if(is.matrix(dist.phylo)){
      D <- dist.phylo[tips, tips]
    } else {
      stop("dist.phylo is not valid")
    }
  }
  if(dist.phylo == "Abouheif"){
    D.max <- max(D[D != max(D)] )
  } else {
    D.max <- max(D)
  }
  
  if(is.null(sigma)){
    sigma <- mean(colMeans(D)/(2 * 1.96))
  }
  
  dist.points <- seq(0, D.max, length.out = n.points)
  res <- matrix(NA, nrow = n.points, ncol = 4)
  for(i in seq(1, n.points)){
    res[i, 1] <- dist.points[i]
    Wi <- phyloWeights(phy, dist.phylo = dist.phylo, method = "lag-norm", mu = dist.points[i], sigma = sigma)
    Wi <- Wi[tips, tips]
    if(n.traits < 2){
      X <- as.vector(t(scale(X)))
      bi <- boot::boot(X, function(x, z) moranTest(xr = x[z], Wr = prop.table(Wi[z, z], 1), reps = 0)$Moran.I, R = ci.bs)
      res[i, 2:3] <- boot::boot.ci(bi, type = "norm", conf = ci.conf)$norm[2:3]
      res[i, 4] <- bi$t0
    } else {
      X <- as.matrix(scale(X))
      bi <- boot::boot(X, function(x, z) mantelStat(xr = x[z, ], Wr = prop.table(Wi[z, z], 1)), R = ci.bs)
      res[i, 2:3] <- boot::boot.ci(bi, type = "norm", conf = ci.conf)$norm[2:3]
      res[i, 4] <- bi$t0
    }
  }
  res[, 2:3][res[, 2:3] > 1] <- 1
  res[, 2:3][res[, 2:3] < -1] <- -1
  
  pcr <- list(res = res, trait = trait, dist.phylo = dist.phylo, sigma = sigma,
              ci.bs = ci.bs, ci.conf = ci.conf, n = n.tips)
  class(pcr) <- "phylocorrelogram"
  return(pcr)  
}


#' Plot a phylogenetic correlogram
#' 
#' This function plots phylogenetic correlograms produced by \code{\link{phyloCorrelogram}}
#' 
#' @param x a \code{phylocorrelogram} object.
#' @param show.ci a logical indicating whether to plot the confidence interval envelop (default \code{TRUE}).
#' @param show.h0 a logical indicating whether to plot the the line showing the expected
#' value of Moran's I under the null hypothesis of no phylogenetic autocorrelation (default \code{TRUE}).
#' @param show.test a logical indicating whether to plot indicator of significance (default \code{TRUE}).
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param main a main title for the plot
#' @param ... other graphical parameters passed to the \code{plot} function.
#' 
#' @examples
#' \dontrun{
#' data(navic)
#' pc <- phyloCorrelogram(navic)
#' plot(pc)
#' }
#' @method plot phylocorrelogram
#' @export
plot.phylocorrelogram <- function(x, show.ci = TRUE, show.h0 = TRUE, show.test =TRUE,
                                  xlab = "Phylogenetic distance", ylab = "Correlation",
                                  main = "Phylogenetic correlogram", ...){
  
  if(class(x) != "phylocorrelogram"){
    stop("x must be an object of class 'phylocorrelogram'.")
  }
  
  h0 <- -1/(x$n-1)
  
  plot(NA, type = "n",
       xlim = range(x$res[, 1]),
       ylim = range(x$res[, 2:3]),
       xlab = xlab, ylab = ylab, main = main, ...)
  #abline(h = 0, col = "grey30")
  if(show.h0){
    abline(h = h0, col = "grey30")
  }
  if(show.ci){
    lines(x = x$res[, 1], y = x$res[, 2], lty = "dashed")
    lines(x = x$res[, 1], y = x$res[, 3], lty = "dashed")
  }
  lines(x = x$res[, 1], y = x$res[, 4], lwd = 2)
  
  if(show.test){
    corpos <- x$res[, 2] > h0 & x$res[, 3] > h0
    corneg <- x$res[, 2] < h0 & x$res[, 3] < h0
    n.res <- dim(x$res)[1]
    bar.col <- rep("black", length.out = n.res)
    bar.col[corpos] <- "red"
    bar.col[corneg] <- "blue"
    segments(x0 = x$res[, 1][-n.res],
             y0 = rep(par("usr")[3], length.out = n.res-1),
             x1 = x$res[, 1][-1],
             y1 = rep(par("usr")[3], length.out = n.res-1),
             col = bar.col, lwd = 10, lend = 1)
  }
  
}

