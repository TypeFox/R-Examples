#' Select markers in a region around QTL
#' 
#' Outputs a list of markers in regions around QTL. Region is defined by the window parameter of x cM to either side of the QTL positions.
#' @export
#' @param qtlpos Vector of QTL positions 
#' @param qtlchr Vector of QTL chromosomes (one for each position)
#' @param map Linkage map to determine where QTL are relative to other markers
#' @param window Number of cM to each side of QTL in which to find nearby markers
#' @param qtlnam Optional vector of names for the QTL
#' @return Returns a map selecting out regions +- window around the QTL positions.
#' @seealso \code{\link[mpMap]{plotlink.map}}
#' @examples 
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl")
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
#' qmap <- qtlmap(summary(mpq.dat)[,2], summary(mpq.dat)[,1], mpq.dat$map)
#' plotlink.map(qmap)


qtlmap <- function(qtlpos, qtlchr, map, window=10, qtlnam) 
{
  ## function to a) plot map in regions around qtl
  ## b) output a list of markers in regions around qtl
  ## region is defined by +- window to either side of qtl positions

  if (length(qtlpos)!=length(qtlchr)) stop("Must input the same number of QTL positions and chromosomes\n")

  ## qtlchr can either be input as a list of names which match the names of map
  ## or as a numeric list indicating the indices of chr in map

  if (missing(qtlnam)) qtlnam <- paste("QTL", 1:length(qtlpos), sep="")
  if (is.numeric(qtlchr)) 
	if (length(setdiff(1:length(map), qtlchr))>0) stop("QTL chromosomes are not contained in map\n") else qtlchr <- names(map)[qtlchr]
  qtlchr <- as.character(qtlchr)

  output <- list()
  ## output is map around each qtl
  for (i in 1:length(qtlpos))
  {
	mapi <- map[[qtlchr[i]]]
	output[[i]] <- mapi[which(abs(mapi-qtlpos[i])<window)]
	## add in the QTL position?
	output[[i]] <- c(qtlpos[i], output[[i]])
	names(output[[i]])[1] <- qtlnam[i]
	output[[i]] <- sort(output[[i]])

	attr(output[[i]], "qtlpos") <- qtlpos[i]
	attr(output[[i]], "qtlchr") <- qtlchr[i]
  }
  names(output) <- qtlchr

  class(output) <- "map"
  output
}

