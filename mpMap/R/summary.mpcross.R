#' Summary of mpcross object
#'
#' Summarizes mpcross object in terms of number of markers, lines, type of markers, and quality of markers
#' @S3method summary mpcross
#' @method summary mpcross
#' @param object Object of class \code{mpcross}
#' @param ... Additional arguments
#' @return Printed summary including - markers which had to be removed due to monomorphic or missing founders; numbers of biallelic/multiallelic markers; percent of markers with missing data; percent of markers with high segregation distortion.
#' @seealso \code{\link[mpMap]{mpcross}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' summary(sim.dat)

summary.mpcross <-
function(object, ...)
{
  dat <- clean(object)

  drop1 <- dat$drop1 
  drop2 <- dat$drop2

  fdr.alleles <- dat$alleles
  pctmiss <- dat$missing

  segpval <- dat$seg[,3]

  cat("-------------------------------------------------------\n")
  cat("Summary of mpcross object\n")
  cat("-------------------------------------------------------\n")
  cat(length(drop2), " markers were removed with missing values in founders\n")
  cat(length(drop1), " markers were removed with non-polymorphic founder genotypes\n")
  cat("-------------------------------------------------------\n")
  cat(length(which(fdr.alleles==2)), " markers were biallelic.\n")
  cat(length(which(fdr.alleles>2) ), " markers were multiallelic.\n")
  cat("-------------------------------------------------------\n")
  cat(length(which(pctmiss>.05)), " markers had >5% missing data.\n")
  cat(length(which(pctmiss>.10)), " markers had >10% missing data.\n")
  cat(length(which(pctmiss>.20)), " markers had >20% missing data.\n")
  cat("-------------------------------------------------------\n")
  cat(length(which(segpval<1e-5)), " markers had <1e-5 p-value for segregation distortion\n")
  cat(length(which(segpval<1e-10)), " markers had <1e-10 p-value for segregation distortion\n")
  cat(length(which(segpval<1e-15)), " markers had <1e-15 p-value for segregation distortion\n")

  
}

