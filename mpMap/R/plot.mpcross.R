#' Plot summary of mpcross object
#'
#' Plots summary of phenotypes and genetic map for mpcross object. If calculated, plots a heatmap of recombination fraction estimates and transformed LOD scores with legend.
#' @S3method plot mpcross
#' @method plot mpcross
#' @param x Object of class \code{mpcross}
#' @param chr Selected chromosomes. Default is all
#' @param ... Additional arguments to plot functions
#' @return Phenotype distributions are plotted as histograms or barplots depending on the variable type. 
#' 
#' Genetic maps are plotted using a version of \code{plotlink.map} for all genetic maps included in the object (may be both a simulated version and an estimated version). 
#' 
#' If recombination fractions have been estimated, a heatmap with recombination fraction estimates below the diagonal and scaled LOD scores above the diagonal is also plotted. LOD scores are transformed to 2^(-(LOD/4+1)) in order to be on the same scale as the theta values.  
#' @seealso \code{\link[mpMap]{mpestrf}}, \code{\link[mpMap]{mpcross.object}}, \code{\link[Heatplus]{heatmap_2}}, \code{\link[mpMap]{plotlink.map}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' plot(sim.dat)

plot.mpcross <- 
function(x, chr, ...) 
{
  par(ask=TRUE)

  noplots <- TRUE
  ### divided up into sections. Some can't be run if parts of object are missing
  if (!is.null(x$map)){
    nmap <- length(grep("map", names(x)))
 
  if (missing(chr)) chr <- c(1:length(x$map))

  ################# This section is the plotted genetic map
  for (i in grep("map", names(x)))
	plotlink.map(x[[i]], marker.names=FALSE)
  noplots <- FALSE
  }

  if (!is.null(x$pheno)) 
  for (i in 1:ncol(as.matrix(x$pheno)))
  {
	x$pheno <- as.data.frame(x$pheno)
	if (is.numeric(x$pheno[,i]))
	  hist(x$pheno[,i], col="violetred", xlab=paste("phe", i), main=colnames(x$pheno)[i]) else if (is.factor(x$pheno[,i]))
	plot(x$pheno[,i], col="royalblue", xlab=paste("phe", i), main=colnames(x$pheno)[i])
	noplots <- FALSE
  }

  ################# This section is for recombination fractions
  if (!is.null(x$rf)){
    require(gdata)
    require(Heatplus)

    markers <- vector()

    for (i in setdiff(1:length(x$map), chr))
      markers <- c(markers, names(x$map[[i]])) 
 
    mpred <- x
    if (length(markers)>0)
  	mpred <- subset(x, setdiff(colnames(x$finals), markers))

    nmrk <- ncol(mpred$founders)
    mat <- matrix(nrow=nmrk, ncol=nmrk)

    upperTriangle(mat) <- 2^(-upperTriangle(mpred$rf$lod)/4-1)
    lowerTriangle(mat) <- lowerTriangle(mpred$rf$theta)
  
    heatmap_2(t(mat), Rowv=NA, Colv=NA, scale="none", legend=1)
    noplots <- FALSE
  }
  if (noplots)
	cat("Object does not contain any plottable components")

}
