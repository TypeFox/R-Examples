#' Detect a second QTL in a QTL profile from (composite) interval mapping
#'
#' Given the output from an initial scan of chromosomes with significant genetic
#' variation, locates the second peak in a QTL profile. 
#' @export
#' @param mpqtl Object of class \code{mpqtl}
#' @param window cM on each side of initial QTL from which to exclude second peak
#' @param drop Drop from original peak required in order to detect a second peak which is above the significance threshold
#' @param dwindow Window over which to smooth p-values - default is five markers
#' @return The original input object with additional entries for newly detected QTL. 
#' @seealso \code{\link[mpMap]{mpIM}}, \code{\link[mpMap]{plot.mpqtl}}, \code{\link[mpMap]{summary.mpqtl}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
#' mpq2 <- findqtl2(mpq.dat, drop=2)
#' plot(mpq2)
#' summary(mpq2)

findqtl2 <- function(mpqtl, window, drop, dwindow=5)
{
  require(VPdtw) 
  output <- mpqtl 
  nchr <- length(mpqtl$QTLresults$qtl)
  support <- supportinterval(mpqtl)$support

  if (missing(window)) window <- apply(support, 2, diff) else if (length(window)==1) window <- rep(window, nchr)

  summ <- summary(mpqtl)  
  chr <- as.character(summ$Chr)
  pkpos <- summ$Pos
  summ$pvalue[summ$pvalue==0] <- 1e-50
  sig <- -log10(summ$pvalue)
  fmap <- attr(mpqtl$prob, "map")
  pkind <- unlist(lapply(1:nchr, function(x) return(which.min(abs(pkpos[x]-fmap[[chr[x]]])))))
  fndrfx <- mpqtl$QTLresults$fndrfx
  se <- mpqtl$QTLresults$se
  pval <- mpqtl$QTLresults$pval
  thr2 <- -log10(attr(mpqtl$QTLresults$qtl, "threshold"))
   

  if (missing(drop)) thr <- rep(thr2, nchr) else thr <- sig-drop

  for (i in 1:nchr)
  {
    cmap <- fmap[[chr[i]]]
	pval[[chr[i]]][pval[[chr[i]]]==0] <- 1e-50
	cpv <- -log10(pval[[chr[i]]])
	cpv2 <- cpv
	cpv[(cmap > pkpos[i]-window[i]) & (cmap < pkpos[i]+window[i])] <- 0 
	stop <- 0

	tmp1 <- dilation(cpv2, dwindow)
	tmp2 <- cpv2==tmp1

	tmp3 <- cbind(cpv2[tmp2], which(tmp2))
	if (sum(tmp2)>1) {
	ind2 <- tmp3[order(tmp3[,1], decreasing=TRUE)[2],2]
	if (cpv2[ind2] > thr2 & min(cpv2[ind2:pkind[i]]) < thr[i])
 	  {
	    k <- attr(output$QTLresults$qtl[[chr[i]]], "index")
  	    output$QTLresults$qtl[[chr[i]]] <- rbind(output$QTLresults$qtl[[chr[i]]], c(cmap[ind2], fndrfx[[chr[i]]][,ind2], se[[chr[i]]][,ind2]))
	    attr(output$QTLresults$qtl[[chr[i]]], "index") <- c(k, ind2)
	  }  
        }
  }
  attr(output$QTLresults$qtl, "nqtl") <- sum(unlist(lapply(output$QTLresults$qtl, nrow)))	  

  output
}
