#' Standardization of landmark configuration by a reference landmark
#'
#' This function provides options for reflecting GPA coordinates around x=0, swapping of GPA xy-coordinates, 
#' and rotating all GPA coordinates in such a way that landmark 7 is standardized at x=0.  
#' @param x an array containing landmark coordinate data of anchor from the specimens of interest
#' @param reflect if TRUE, x-coordinates are reflected around the x=0 axis
#' @param swap if TRUE, the x and y-coordinates are swapped
#' @param sgn a numeric vector; two choices 1 and -1; defaults to c(1,1)
#' @return an array of landmark coordinates
#' @details The values in the \code{sgn} vector should be tuned to obtain the desired orientation of landmark configuration.
#' This function processes the output from \code{procrustesFit}, which requires the \code{geomorph} package 
#' (Version 3.0.0; Adams & Otarola-Castillo, 2013). 
#' @seealso \code{\link{procrustesFit}}, \code{\link{plotLM}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Adams DC, Otarola-Castillo E. (2013). geomorph: an R package for the collection and analysis of geometric morphometric shape data. 
#' Methods in Ecology and Evolution 4:393-399.
#'
#' Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' library(geomorph)
#'
#' data(ligophorus_tpsdata)
#' #A data processing step to parse out the orientation of landmarks 
#' #from samples of L.parvicopulatrix
#'
#' O <- matrix(0, length(ligophorus_tpsdata$parvicopulatrix), 4)
#' for(w in 1:length(ligophorus_tpsdata$parvicopulatrix)){
#'	result <- mapply(function(k)
#'	anglePolygon(matrix2list(ligophorus_tpsdata$parvicopulatrix[[w]][(11*(k-1)+1):(11*k),]),
#'	degree=TRUE), k=1:4)
#'
#'	result_angle <- mapply(function(k) list(result[[2*k-1]]), k=1:4)
#'	result_orientation <- mapply(function(k) list(result[[2*k]]), k=1:4)
#'	names(result_angle) <- names(result_orientation) <- c("VR","VL","DR","DL")
#'	O[w,] <- unlist(result_orientation)
#' }
#'
#' mdir <- apply(O, 2, function(k) which(k == "m") )
#' pdir <- apply(O, 2, function(k) which(k == "p") )
#'
#' e <- 1 #Ventral right anchor
#' result <- procrustesFit(ligophorus_tpsdata$parvicopulatrix, e, 
#' list(mdir[[e]], pdir[[e]]), PrinAxes=TRUE, showplot=TRUE)
#' #Standardize the x-coordinate of Landmark 7 by rotating the x-coordinate 
#' #of its mean GPA xy-coordinate to x=0.
#' coordinates <- stdLM(result$coords, reflect=FALSE, swap=TRUE, sgn=c(1,-1))
#'
#' plotLM(coordinates, "VR", pointscale=0.8,axispointscale=0.8,
#' meansize=1.2,polygon.outline=TRUE,c(-.6,.6),c(-.6,.6) )
#'

stdLM <- function(x,reflect=FALSE,swap=FALSE,sgn=c(1,1)){

VV <- x
PP <- VV

if(reflect == TRUE){
	for(i in 1:dim(VV)[3]) {
		VV[,,i][,1] <- VV[,,i][,1] * -1
	}
	PP <- VV
}

if(swap == TRUE){
	for(i in 1:dim(VV)[3]){
		PP[,,i][,1] <- sgn[1]*VV[,,i][,2]
		PP[,,i][,2] <- sgn[2]*VV[,,i][,1]
	}
}

	angle.displace <- function(y, ref){
	refland <- ref
	radius0 <- sqrt( sum ((y[refland,])^2 ) )
	theta0 <- asin( abs(y[refland,1]) / radius0 )
	return(theta0)
	}

	rotateRef <- function(x, theta0){

	xsign <- t(apply(x, 1, sign))

	radius <- apply(x[,], 1, function(k) sqrt ( sum(k^2) ) )
	P <- apply(xsign,1,prod)

	#radius is the Euclidean distance between a landmark and the origin
	omega <- acos( abs(x[,1]) / radius )
	xnew <- radius * cos(omega + theta0) * (P==1) +  radius * cos(abs(omega - theta0) ) * (P==-1)
	ynew <- radius * sin(omega + theta0) * (P==1) +  radius * sin(abs(omega - theta0) ) * (P==-1)

	return(cbind(xnew, ynew) * xsign)

	}

G <- PP
GG <- G[,,1]

for(i in 2:dim(G)[3]){
	GG <- GG + G[,,i]
}

M <- GG / dim(G)[3]

refland <- 7

refsgn <- sign ( M[refland,1] )
theta <- refsgn * (-1) * angle.displace(M,refland) 

Glist <- G
for(i in 1:dim(G)[3]){
	Glist[,,i] <- rotateRef(G[,,i],theta)
}

return(Glist)

}