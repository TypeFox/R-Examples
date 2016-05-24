#' Summary of mpprob object
#' 
#' Summarizes details about underlying mpcross object as well as descriptive statistics about estimated founder haplotypes
#' @S3method summary mpprob
#' @method summary mpprob
#' @param object Object of class \code{mpprob}
#' @param ... Additional arguments
#' @return Output to screen of percentage of each chromosome inherited from founders, average number of recombinations per chromosome and genomewide, and number of finals/founders/chromosomes/markers per chromosome. 
#' @seealso \code{\link[mpMap]{plot.mpprob}}, \code{\link[mpMap]{mpprob}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl")
#' summary(mpp.dat)

summary.mpprob <- function(object, ...)
{
  n.founders <- nrow(object$founders)
  ## summary of how it was computed - threshold/step size
  cat("This is a summary of an mpprob object.\n")
  cat("Probabilities were computed at a step size of ", attr(object$prob, "step"), "\n")
  cat("Founders were estimated using a threshold of ", attr(object$estfnd, "threshold"), "\n")
  cat("----------------------------------------------------------\n\n")
  cat("Percent of each chromosome with each founder ancestry (NA indicates unknown)\n")

  ## counting the % of the genome with each type of ancestry
  cts1 <- lapply(object$estfnd, function(x) {
	 	z <- factor(as.vector(x), levels=1:n.founders)
		return(round(table(z, useNA="always")/prod(dim(x))*100, 2))})
  cts <- do.call("cbind", cts1)

  ### count the number of recombinants per line, per chromosome, then sum up, and average per 100 cM
  nrec1 <- lapply(object$estfnd, function(x) return(apply(x, 1, function(y) return(sum(diff(y[!is.na(y)])!=0)))))
  nrec <- do.call("cbind", nrec1)
  chrlength <- unlist(lapply(object$map, max))/100
  genlength <- sum(chrlength)
  avgrec <- round(mean(apply(nrec, 1, sum))/genlength, 2)
  avgrecchr <- round(apply(nrec, 2, mean)/chrlength, 2)
  names(avgrecchr) <- names(object$map)

  colnames(cts) <- names(object$map)

  print(cts)
  
  cat("\n----------------------------------------------------------\n")
  cat("Average number of recombinations per chr (per M):\n")
  print(avgrecchr)
  cat("\nTotal genome length (in M): ", genlength, "\n")
  cat("Average number of recombinations (per M): ", avgrec, "\n")

  cat("----------------------------------------------------------\n")
  cat("This object is based off of an mpcross with:\n")
  cat("No. of finals:       ", nrow(object$finals), "\n")
  cat("No. of founders:     ", nrow(object$founders), "\n")
  cat("No. of chromosomes:  ", length(object$map), "\n")
  cat("Markers per chr:     ", unlist(lapply(object$map, length)), "\n")
}

