#' Mass Extraction of Generalized Procrustes Analysis Coordinates from Anchors
#'
#' Given a list of landmark coordinate data, this function performs Generalized Procrustes Analysis (GPA) and extracts the GPA-coordinates. 
#  from all ventral and dorsal anchors   
#' @param x a list containing landmark coordinate data of anchors from the specimens of interest
#' @param e a constant specifying the anchor of interest: ventral right(1), ventral left(2), dorsal right(3), dorsal left(4)
#' @param makeplot if TRUE, returns a scatter plot of the GPA-coordinates
#' @param reflect logical; if TRUE, x-coordinates are reflected around the x=0 axis
#' @param swap logical; if TRUE, the x and y-coordinates are swapped
#' @param axispointscale a numeric constant for controlling the font size of numeric values on the xy axes
#' @param sgn a numeric vector; two choices 1 and -1; defaults to c(1,-1)
#' @return an array containing GPA-coordinates of the specimens of interest
#' @details This function is essentially a wrapper for \code{procrustesFit} and \code{stdLM} to ease extraction of GPA-coordinates
#' from list data. Both require the \code{geomorph} package (Version 3.0.0).
#' @seealso \code{\link{procrustesFit}}, \code{\link{stdLM}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Adams DC, Otarola-Castillo E. (2013). geomorph: an R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution 4:393-399.
#' @examples
#' data(ligophorus_tpsdata)
#'
#' vright <- procrustesFit.2(ligophorus_tpsdata$johorensis, 1,
#' makeplot=TRUE, reflect=FALSE, swap=TRUE, sgn=c(-1,1))
#'
#' vleft <- procrustesFit.2(ligophorus_tpsdata$johorensis, 2, 
#' makeplot=TRUE, reflect=TRUE, swap=TRUE, sgn=c(-1,1))
#'
#' va <- (vright+vleft)/2
#'
#' plotLM(va, "VA", pointscale=0.8, meansize=1.2, polygon.outline=TRUE, 
#' axispointscale=0.8, c(-.6,.6),c(-.6,.6))

procrustesFit.2 <- function(x,e,makeplot=FALSE, reflect=FALSE, swap=TRUE, axispointscale=0.8, sgn=c(1,-1)){

O <- matrix(0,length(x),4)
type <- c("VR","VL","DR","DL")

#A xa processing step to parse out the orientation of landmarks for each specimen
for(j in 1:length(x)){
	result <- mapply(function(k)
	anglePolygon(matrix2list(x[[j]][(11*(k-1)+1):(11*k),]),degree=TRUE),
	k=1:4)

	result_angle <- mapply(function(k) list(result[[2*k-1]]), k=1:4)
	result_orientation <- mapply(function(k) list(result[[2*k]]), k=1:4)
	names(result_angle) <- names(result_orientation) <- type
	O[j,] <- unlist(result_orientation)
}

mdir <- apply(O, 2, function(k) which(k == "m") )
pdir <- apply(O, 2, function(k) which(k == "p") )

result <- procrustesFit(x, e, list(mdir[[e]], pdir[[e]]), PrinAxes=TRUE, showplot=FALSE)
coordinates <- stdLM(result$coords, reflect=reflect, swap=swap, sgn=sgn)

if(makeplot==TRUE){
plotLM(coordinates, type[e], pointscale=0.8, axispointscale=axispointscale, meansize=1.2,polygon.outline=TRUE,c(-.6,.6),c(-.6,.6) )
}

return(coordinates)
}
