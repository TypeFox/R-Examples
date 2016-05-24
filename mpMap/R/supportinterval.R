#' Calculate support interval for detected QTL
#'
#' Calculates support interval for QTL based on Wald profile and QTL position
#' @export
#' @param x Object of class \code{mpqtl}
#' @param chr Selected chromosomes
#' @param lodsupport Size of support interval; default is 1 LOD
#' @details Computes the x-LOD support interval as the region surrounding a QTL peak in  which the Wald profile exceeds the equivalent of x LOD less than the peak value. 
#' @return A list with two components: a matrix containing lower and upper bounds for the support intervals for each QTL, and the positions of QTL on each chromosome. 
#' @seealso \code{\link[mpMap]{plot.mpqtl}}
#' @examples
#' map <- sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
#' si <- supportinterval(mpq.dat)

supportinterval <- function(x, chr, lodsupport=1)
{
  map <- attr(x$prob, "map")
  if (missing(chr)) chr <- 1:length(map)
  if (is.factor(chr)) chr <- as.character(chr)
  if (is.numeric(chr)) chr <- names(map)[chr]

  ## set up as a scanone output
  waldsc <- matrix(nrow=length(unlist(x$QTLresults$wald)), ncol=3)
  waldsc[,1] <- rep(names(x$QTLresults$wald), unlist(lapply(x$QTLresults$wald, length)))
  waldsc[,2] <- unlist(map)
  waldsc[,3] <- unlist(x$QTLresults$wald)

  waldsc <- as.data.frame(waldsc)
  waldsc[,2] <- as.numeric(as.character(waldsc[,2])) 
  waldsc[,3] <- as.numeric(as.character(waldsc[,3])) 

  names(waldsc) <- c("chr", "pos", "lod")
  class(waldsc) <- c("scanone", "data.frame")

  map2 <- map[chr]
  qtlpos <- vector()
  chrpos <- c(0,cumsum(unlist(lapply(map2, max))+25))
  for (i in match(names(x$QTLresults$qtl), names(map2)))
	qtlpos <- c(qtlpos, x$QTLresults$qtl[[names(map2)[i]]][,1]+chrpos[i])

### instead of plotting a line to indicate the QTL; plot grey region between
  ### two lines to indicate support interval

  ### determine appropriate LOD threshold by converting relative df
  ### if number of QTL detected > 0
  output <- list()
  if (length(x$QTLresults$qtl)>0) {
  qtlmat <- do.call("rbind", x$QTLresults$qtl)
  rownames(qtlmat) <- rep(names(x$QTLresults$qtl), unlist(lapply(x$QTLresults$qtl, nrow)))

  sil <- vector()
  siu <- vector()
  thr <- qchisq(pchisq(lodsupport*3.84,1), 3)
  chrnam <- vector()
  for (i in 1:length(qtlpos)) {
	 m <- match(rownames(qtlmat)[i], names(map2))
	 chrstart <- nrow(waldsc) - max((nrow(waldsc)-(1:nrow(waldsc)))*(waldsc[,1]==rownames(qtlmat)[i]))
	 chrend <- max((1:nrow(waldsc))*(waldsc[,1]==rownames(qtlmat)[i]))
	 waldsc2 <- waldsc[chrstart:chrend,]
	 peakind <- which(abs(waldsc2[,2]-qtlmat[i,1])<.001)
	if (length(peakind)>1) peakind <- peakind[1]
#	 sil[i] <- waldsc[max(c(chrstart, (chrstart:peakind)*(waldsc[chrstart:peakind, 3] < (waldsc[peakind,3]-thr)))),2]
#	 siu[i] <- waldsc[chrend-max((chrend-peakind:chrend)*(waldsc[peakind:chrend,3] < (waldsc[peakind,3]-thr))),2]
	sil[i] <- waldsc2[max(c(1, (1:peakind)*(waldsc2[1:peakind,3] < (waldsc2[peakind,3]-thr)))), 2]
	siu[i] <- waldsc2[nrow(waldsc2)-max((nrow(waldsc2)-peakind:nrow(waldsc2))*(waldsc2[peakind:nrow(waldsc2), 3] < (waldsc2[peakind,3]-thr))),2]

#	sil[i] <- min(waldsc2[(waldsc2[,3] > (waldsc2[peakind,3]-thr)),2])
#	siu[i] <- max(waldsc2[(waldsc2[,3] > (waldsc2[peakind,3]-thr)),2])
  }

  output$support <- rbind(sil, siu)
  rownames(output$support) <- c("Lower", "Upper")
  colnames(output$support) <- chrnam
  
  chrnam <- rownames(qtlmat)
  names(qtlpos) <- chrnam
  output$qtlpos <- qtlpos
  }
  output
}
