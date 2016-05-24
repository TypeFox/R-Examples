#' Estimate pairwise recombination fractions between markers
#'
#' Estimates pairwise recombination fractions by maximizing the likelihood for a multi-parent cross over a grid of possible values. Theta values and corresponding LOD scores are returned for each pair of markers in the object.
#' @export
#' @useDynLib mpMap
#' @param object Object of class \code{mpcross}
#' @param r Grid of potential recombination values. If missing the function will maximize over (0, .005, .01, .015, ... , .095, .1, .11, .12, ... .49, .5). 
#' @param grid Flag for whether to output the entire grid of likelihoods at each potential recombination value
#' @return Returned object is of the class 'mpcross' with the additional component \code{rf}. If n.mrk is the number of markers genotypes, this is a list with components:
#' \item{rf$theta}{ n.mrk x n.mrk matrix of estimated recombination fractions between each pair of loci}
#' \item{rf$lod}{ n.mrk x n.mrk matrix of LOD scores at the estimated recombination values}
#' \item{rf$lkhd}{ n.mrk x n.mrk matrix of likelihood values at the estimated recombination values}
#' @seealso \code{\link[mpMap]{mpcross}}, \code{\link[mpMap]{plot.mpcross}}
#' @examples
#' map <- sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, qtl=matrix(data=c(1, 50, .4, 0, 0, 0), nrow=1, ncol=6, byrow=TRUE), seed=1)
#' dat.rf <- mpestrf(sim.dat)
#' plot(dat.rf)

mpestrf <-
function(object, r, grid=FALSE)
{ 
  if (!inherits(object, "mpcross")) stop("Object must be of class mpcross")

  require(gdata)
 
  if (missing(object))	
	stop("Missing a required argument for this function")

  if (missing(r)) r <- c(0:20/200, 11:50/100)

  output <- object 
  founder.g <- object$founders
  final.g <- cbind(object$id, object$finals)
  pedigree <- object$pedigree

  n.founders <- nrow(founder.g)
  n.loci <- ncol(founder.g)
  n.finals <- nrow(final.g)

  # pairs of loci
  if (n.loci>1) {
    pairsij <- t(cbind(combn(1:n.loci, 2), matrix(data=rep(1:n.loci,2), nrow=2, ncol=n.loci, byrow=TRUE)))
  pairsij <- pairsij[order(pairsij[,1], pairsij[,2]), c(2:1)]
  } else pairsij <- matrix(data=rep(1:n.loci,each=2), nrow=n.loci, ncol=2, byrow=TRUE)

  n.pairs <- nrow(pairsij)

  rpairs <- CR_estrf(final.g, founder.g, pedigree, pairsij, r)

  if (grid)  output$lkhdgrid <- cbind(pairsij, rpairs)

  output$rf <- list()
  output$rf$theta <- matrix(nrow=n.loci, ncol=n.loci)
  output$rf$lod <- matrix(nrow=n.loci, ncol=n.loci)
  output$rf$lkhd <- matrix(nrow=n.loci, ncol=n.loci)

  maxlkhd <- apply(rpairs, 1, max)
  minlkhd <- apply(rpairs, 1, min)

  same <- ((maxlkhd-minlkhd)<.000001)
  tm <- r[apply(rpairs, 1, which.max)]
  tm[which(same==TRUE)] <- NA
  lowerTriangle(output$rf$theta, diag=TRUE) <- tm

  theta.5 <- which(r==0.5)
  lm <- apply(rpairs, 1, max)
  lowerTriangle(output$rf$lkhd, diag=TRUE) <- lm 
  lm <- apply(rpairs-rpairs[,theta.5], 1, max)
  lm[which(same==TRUE)] <- NA
  lowerTriangle(output$rf$lod, diag=TRUE) <- lm

  if (grid) output$lodgrid <- cbind(pairsij, rpairs-rpairs[,theta.5])

  upperTriangle(output$rf$lod) <- upperTriangle(t(output$rf$lod))
  upperTriangle(output$rf$theta) <- upperTriangle(t(output$rf$theta))
  upperTriangle(output$rf$lkhd) <- upperTriangle(t(output$rf$theta))

  rownames(output$rf$lod) <- rownames(output$rf$theta) <- rownames(output$rf$lkhd) <- colnames(output$rf$lod) <- colnames(output$rf$theta) <- colnames(output$rf$lkhd) <- colnames(object$finals)

  return(output)
}

