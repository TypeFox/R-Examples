#' Internal angles of a polygon
#'
#' This function computes the internal angles of a polygon, given the coordinates of its vertices.
#' @param point a list containing coordinates of the polygon vertices in either clockwise or anti-clockwise direction
#' @param degree if TRUE, returns the internal angles in degrees instead of radians
#' @return A list with two components:
#' \item{angle}{a matrix; the rows represents the vertices and the columns contain the latters' geometrical information:  
#' first column for x coordinate, second column for y coordinate, and third column for the assocaited internal angle (in radians)}  
#' \item{orientation}{a character indicating the direction of ordering the vertices: \code{m} for anti-clockwise and \code{p} for clockwise} 
#' @details The sum of all the internal angles of a polygon with \emph{n}-vertices must be equal to the product of 
#' \emph{n-2} with 180 (degrees) or pi (radians).
#' This function is useful for detecting tps data files that contain errors (e.g. wrong sequence of digitizing landmarks, missing landmarks)
#' so that corrective steps can be taken.
#' @author Thian Liang Cheow \email{Tl2cheow@@yahoo.com}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' #internal angles of a right-angled triangle in degrees
#' anglePolygon(list(c(0,1),c(0,0),c(1,0)),degree=TRUE)
#'
#' data(ligophorus_tpsdata)
#' 
#' #polygonal approximation of anchor shape
#' #the right ventral anchor
#' anchorpolygon <- anglePolygon(matrix2list(ligophorus_tpsdata$bantingensis[[1]][1:11,]), degree=TRUE)
#'
#' #there are 11 landmarks, do the sum of internal angles should be (11-2)x180 = 1620
#' sum(anchorpolygon$angle[,3])
#' 
#' #does this make sense?
#' polyVis(1, havelist=TRUE, listdata=ligophorus_tpsdata$bantingensis)
#'

# Find internal angles of a given polygon, given a list of point in the form list(c(x1,y1), c(x2,y2), ...))
# REQUIREMENT: The list of point must be given in counter clockwise direction. 
# if direction is clockwise, the external angles will be returned.
# Procedures
# 1. For each point, p[i], use it as the center of polar coordinate to find the angle of its two 
# neighbouring point, p[i-1] and p[i+1]. Assume that a1 and a2 are the result. 
# Take p[0] = p[n] and p[n + 1] = p[0] where n is the length of the list.
# 2. Find the difference between a1 and a2, a1 - a2. If the result is negative add 2pi.

anglePolygon <- function(point, degree=FALSE) {

# Get the angle of a given point, range 0 and 2pi
findAngle <- function(p) {
	return (atan2(p[2], p[1]) + (p[2]<0)*2*pi)
}
	count <- length(point)
	result = list()
	# For each point in the list, find the angles of the edges connecting current and neighbouring point.
	# Find the difference of the angles.
	for (i in 1:count) {
		p <- point[[i]]
		if (i == 1) {
			p1 <- point[[count]]
			p2 <- point[[i + 1]]			
		} else if (i == count) {
			p1 <- point[[i - 1]]
			p2 <- point[[1]]			
		} else {
			p1 <- point[[i - 1]]
			p2 <- point[[i + 1]]
		}

		# Need to look into reason of getting negative value
		a <- findAngle(p1 - p) - findAngle(p2 - p)
		a <- a + (a<0)*2*pi
		if(degree==TRUE){
		a <- a/(2*pi) * 360
		}
		result <- append(result, a)		
	}

	output <- mapply(c, point, result, SIMPLIFY=FALSE)	
	final_output <- matrix(unlist(output),ncol=3,byrow=TRUE)
	colnames(final_output) <- c("x","y","angle")
     total_angle <- (degree==TRUE) *round(sum(final_output[, 3]), 0) + (degree==FALSE)*round(sum(final_output[, 3])/2/pi * 360, 0)
	direction <- "m"

	if(degree == TRUE){
	wrong_orientation <- total_angle > 180*(nrow(final_output)-2) 
	}

	else if(degree == FALSE){ 
	wrong_orientation <- total_angle/360 * 2 * pi != pi*(nrow(final_output)-2) 
	}

	if(wrong_orientation == TRUE){
	final_output[,3] <- (degree==TRUE)*(360-final_output[,3]) + (degree==FALSE)*(2*pi-final_output[,3])
	direction <- "p"
	}

	out <- list(final_output, direction)
	names(out) <- c("angle","orientation")
	return(out)
}
