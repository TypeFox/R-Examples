#' Plot summary of founder probabilities and haplotype blocks
#' 
#' Plot the percentage of each chromosome inherited from each founder
#' @S3method plot mpprob
#' @method plot mpprob
#' @param x Object of class \code{mpprob}
#' @param chr Chromosomes to plot. Default is all
#' @param nlines Number of lines to select to show founder ancestry. Default is all
#' @param ... Additional arguments to plot function
#' @return Barplot of the percentage of each founder on each chromosome; individual heatmaps of which chunks of each chromosome are inherited from each founder.#' @seealso \code{\link[mpMap]{mpprob}}, \code{\link[mpMap]{summary.mpprob}}, \code{\link[Heatplus]{heatmap_2}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl")
#' plot(mpp.dat)

plot.mpprob <-
function(x, chr, nlines, ...)
{
  if (missing(chr)) chr <- 1:length(x$map)
  if (missing(nlines)) nlines <- nrow(x$finals)
  n.founders <- nrow(x$founders)

  cts1 <- lapply(x$estfnd, function(y) {
	z <- factor(as.vector(y), levels=1:n.founders)
	return(round(table(z)/prod(dim(y))*100, 2)) })
  cts <- do.call("cbind", cts1)

  ## first plot is the stacked barplot of founder probabilities
  cts1 <- cbind(cts, rep(NA, nrow(cts))) 
  cts2 <- as.matrix(c(rep(100, ncol(cts)), NA))

  barplot(t(cts2), col="white")

  barplot(cts1, col=c("lightblue2", "royalblue1", "darkseagreen1", "seagreen3", "lemonchiffon", "darkgoldenrod1", "indianred1", "firebrick")[1:n.founders], main="Founder %age by Chromosome", xlab="Chromosome", legend=rownames(cts1), add=TRUE)

  par(ask=TRUE)
  sum <- vector(length=nrow(x$finals))
  for (i in chr)
  {
    nrec <- apply(x$estfnd[[i]], 1, function(x) return(sum(diff(x[!is.na(x)])!=0)))
    ## if lines is not missing, need to select out the ones with the most
    ## recombination events
    ord <- cbind(1:length(nrec), nrec)
    ord <- ord[order(ord[,2], decreasing=TRUE), ]    
    if (length(nlines)>1) sel <- nlines else sel <- ord[1:nlines, 1]

    mat <- x$estfnd[[i]][sort(sel),]
    colnames(mat) <- names(attr(x$prob, "map")[[i]])
    rownames(mat) <- rownames(x$finals)[sort(sel)]

    heatmap(t(mat), col=c("lightblue2", "royalblue1", "darkseagreen1", "seagreen3", "lemonchiffon", "darkgoldenrod1", "indianred1", "firebrick")[1:n.founders], main=names(x$map)[[i]], Rowv=NA, Colv=NA, scale="none", xlab="Lines", ylab="Markers")
  }
}

