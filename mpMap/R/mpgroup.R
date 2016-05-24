#' Group markers into linkage groups given 2-pt recombination fraction estimates
#'
#' Groups markers based on estimated pairwise recombination fractions. Linkage groups are built up by adding a marker if it satisfies the criteria for linkage with at least one other marker in the group.
#' @export
#' @param mpcross  an object of class \code{mpcross} which includes the component \code{rf} output by \code{mp.est.rf}. See \code{\link[mpMap]{mpcross.object}} for further details.
#' @param theta Threshold for grouping based on recombination fraction value
#' @param LOD Threshold for grouping based on LOD score value
#' @return The original object, with the added component \code{lg} which is a list with the following components:
#' \item{n.groups}{ The number of linkage groups formed by the function}
#' \item{groups}{ Vector labelling each marker by assigned linkage group. Missing values mean that the marker was linked to more than one group and could not be assigned with confidence}
#' \item{LODthresh}{ The LOD threshold value used to determine linkage}
#' \item{thetathresh}{ The theta threshold value used to determine linkage}
#' \item{order}{ A list with a component for each constructed linkage group which contains the order of the markers within the linkage group}
#' @seealso \code{\link[mpMap]{mpestrf}}
#' @seealso 
#' map <- sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, qtl=matrix(data=c(1, 50, .4, 0, 0, 0), nrow=1, ncol=6, byrow=TRUE), seed=1)
#' dat.rf <- mpestrf(sim.dat)
#' dat.lg <- mpgroup(dat.rf)

mpgroup <-
function(mpcross, theta=.15, LOD=5)
{
  if (missing(mpcross)) 
	stop("Must input mpcross object to create linkage groups")

  if (is.null(mpcross$rf))
	stop("Must calculate recombination fractions prior to grouping loci")

  output <- mpcross

  nmar <- nrow(mpcross$rf$theta)
  groups <- rep(NA, nmar)
  group.count <- 1

  matt <- mpcross$rf$theta
  matl <- mpcross$rf$lod

  # make sure matrix is symmetric
  for (i in 1:nmar)
  for (j in 1:i)
  {
	matt[j,i] <- matt[i,j]
	matl[j,i] <- matl[i,j]
  }
  matt[is.na(matt)] <- .5
  matl[is.na(matl)] <- 0
  diag(matt) <- NA
  diag(matl) <- NA

  for (m in 1:nmar) {
    if (is.na(groups[m])) {
	grouping <- c(m, which(matt[m,]<theta & matl[m,] >LOD))
	if (length(grouping)>1) {
	  flag <- 1
	  while (flag) {
	    flag <- 0
	    for (i in grouping) {
	 	group_parc <- which(matt[i,]<theta & matl[i,]>LOD)
		for (j in group_parc)
	 	if (all(grouping != j)) {
		  grouping <- c(grouping, j)
		  flag <- 1 
		}
	    }
	  }
	  groups[grouping] <- group.count
	  group.count <- group.count+1
	}
    }
  }	

  output$lg <- list(n.groups=max(groups, na.rm=TRUE), groups=groups, LODthresh=LOD, thetathresh=theta)

  output
}

