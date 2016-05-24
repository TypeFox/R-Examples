#'Leaf dispersion of 3D plants
#'
#'@description This function calculates the leaf dispersion for 3D plants, following Duursma et al. (2012).
#'
#'The method is based on the mean distance to k nearest neighbors in 3D.  The
#'function 'leafdispersion' computes this observed mean distance (\code{Ok})
#'for a plant (an object of class \code{plant3d}), as well as for a square box
#'with randomly distributed leaves at the same leaf area density.
#'
#'
#'@param plant An object of class 'plant3d'.
#'@param kneighbors Number of neighbors to be used.
#'@param nreplicate For the random distribution, the number of replicates to
#'simulate.
#'@param crownvol,nleaves Crown volume and number of leaves - optional. If not
#'provided, they are calculated from the 'plant' object.
#'@return A list with the following components: \describe{
#'\item{list("Ok")}{Observed distance to k nearest neighbors}
#'\item{list("Ek_noedge")}{Expected distance to k nearest neighbors, for random
#'distribution; no edge correction.} \item{list("Ek_edge")}{As above, but with
#'an edge correction} \item{list("Ek_edgeSD")}{Standard deviation among
#'replicates of \code{Ek_edge}} \item{list("kneighbors")}{Number of neighbors
#'for distance calculation} \item{list("disp_edge")}{Edge-corrected leaf
#'dispersion (as in Duursma et al. 2012).} \item{list("disp_noedge")}{Non
#'edge-corrected leaf dispersion} }
#'@author Remko Duursma
#'@references Duursma, R.A., D.S. Falster, F. Valladares, F.J. Sterck, R.W.
#'Pearcy, C.H. Lusk, K.M. Sendall, M. Nordenstahl, N.C. Houter, B.J. Atwell, N.
#'Kelly, J.W.G. Kelly, M. Liberloo, D.T. Tissue, B.E. Medlyn and D.S.
#'Ellsworth. 2012. Light interception efficiency explained by two simple
#'variables: a test using a diversity of small- to medium-sized woody plants.
#'New Phytologist. 193:397-408.
#'@keywords misc
#'@examples
#'
#'
#'# Leafdispersion for the Toona plant
#'leafdispersion(toona)
#'
#'@export
leafdispersion <- function(plant, kneighbors=5, nreplicate=10, nleaves=NA, crownvol=NA){
	# Nearest neighbors (n=5) : Ok is observed mean distance to 5 nearest neighbor leaf midpoints,
	# Ek is expected (based on numerical simulations).
	
	if(class(plant) != "plant3d")stop("Need object of class 'plant3d'")
	if(any(sapply(plant$leaves, function(x)any(is.na(x$XYZ)))))return(NA)
		
	if(is.na(nleaves)){
		nleaves <- plant$nleaves
	}

	if(nleaves < min(kneighbors))return(NA)

	if(is.na(crownvol)){
		ch <- crownhull(plant, plotit=FALSE)
		crownvol <- ch$crownvolume * 1E-09  # m3
	}
	
	lad <- nleaves / crownvol  # m-3
  
  # Determined from numerical simulation
	Ek_noedge <- lad^(-1/3) * (0.5 -0.08784436 + 0.21374167 * sqrt(5))
		
	if(length(kneighbors) > 1)
		kneighbors <- kneighbors[1]:min(kneighbors[length(kneighbors)],nleaves-1)
	
	l <- list()
	l$Ok <- Kn(plant,kneighbors)
	
	Ek_edge <- replicate(nreplicate, KE(lad, nleaves, kneighbors))
	
	l$Ek_noedge <- Ek_noedge
	
	if(length(kneighbors) == 1){
		l$Ek_edge <- mean(Ek_edge)
		l$Ek_edgeSD <- sd(Ek_edge)
	} else {
		l$Ek_edge <- apply(Ek_edge,1,mean)
		l$Ek_edgeSD <- apply(Ek_edge,1,sd)
	}
	
	l$kneighbors <- kneighbors
	l$disp_edge <- l$Ok / l$Ek_edge
	l$disp_noedge <- l$Ok / l$Ek_noedge
	return(l)
}

