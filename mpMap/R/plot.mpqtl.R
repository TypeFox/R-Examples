#' Plot output from interval mapping with detected QTL
#'
#' Plot -log10(p-value) or test statistic against cM position for (composite) interval mapping in multi-parent crosses. QTL support intervals are indicated with rectangles surrounding peaks. 
#' @S3method plot mpqtl
#' @method plot mpqtl
#' @param x Object of class \code{mpqtl}
#' @param wald Flag for whether to plot the Wald statistic or -log10(p)
#' @param chr Set of chromosomes to plot
#' @param lodsupport x-LOD support interval plotted in green
#' @param ... Additional arguments to plotting function
#' @return Plots the -log10(p) or Wald statistic for all chromosomes against the total genome in cM. QTL support intervals are indicated with shaded rectangles surrounding peaks 
#' @seealso \code{\link[mpMap]{mpIM}}, \code{\link[mpMap]{summary.mpqtl}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
#' plot(mpq.dat)

plot.mpqtl <- function(x, wald=FALSE, chr, lodsupport=1, ...)
{
  map <- attr(x$prob, "map")
  if (missing(chr)) chr <- 1:length(map)
  if (is.factor(chr)) chr <- as.character(chr)
  if (is.numeric(chr)) chr <- names(map)[chr]

  ## set up as a scanone output
  psc <- data.frame(chr=rep(names(x$QTLresults$wald), unlist(lapply(x$QTLresults$wald, length))), pos=unlist(map), lod=-log10(unlist(x$QTLresults$pval)))
  waldsc <- data.frame(chr=rep(names(x$QTLresults$wald), unlist(lapply(x$QTLresults$wald, length))), pos=unlist(map), lod=unlist(x$QTLresults$wald))

  class(psc) <- class(waldsc) <- c("scanone", "data.frame")

  if (any(psc[,3]==Inf)) { 
	wald <- TRUE
	cat("Some p-values=0; plotting Wald score\n")
    } 

  if (wald==TRUE) qtl:::plot.scanone(waldsc, chr=chr, ylab="Wald") else qtl:::plot.scanone(psc, chr=chr, ylab="-log10(p)")

  map2 <- map[chr]
  qtlpos <- vector()
  chrpos <- c(0,cumsum(unlist(lapply(map2, max))+25))
  for (i in match(names(x$QTLresults$qtl), names(map2)))
	qtlpos <- c(qtlpos, x$QTLresults$qtl[[names(map2)[i]]][,1]+chrpos[i])

  si <- supportinterval(x, chr=chr, lodsupport=lodsupport)$support

  output <- list()
  if (length(x$QTLresults$qtl)>0) {
  qtlmat <- do.call("rbind", x$QTLresults$qtl)
  rownames(qtlmat) <- rep(names(x$QTLresults$qtl), unlist(lapply(x$QTLresults$qtl, nrow)))
  for (i in 1:length(qtlpos)) {
  	chrnam <- vector()
	 m <- match(rownames(qtlmat)[i], names(map2))
	 chrstart <- nrow(waldsc) - max((nrow(waldsc)-(1:nrow(waldsc)))*(waldsc[,1]==rownames(qtlmat)[i]))
	 chrend <- max((1:nrow(waldsc))*(waldsc[,1]==rownames(qtlmat)[i]))

	rect(xleft=si[1,i]+chrpos[m], xright=si[2,i]+chrpos[m], ybottom=-1000, ytop=1000, col="#00990022")
  }

  output$support <- si
  }

  invisible(output)
}
