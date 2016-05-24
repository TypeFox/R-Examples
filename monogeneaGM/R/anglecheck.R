#' Angle swept by a ray 
#'
#' This function computes angle swept by a ray in anti-clockwise direction.
#' The ray starts at the half-line pointing towards the positive x direction,
#' and passes through the origin and a point defined by the user.
#' @param x a vector of length 2 giving the xy-coordinates of the point defined by the user
#' @param radians if FALSE, returns angle in degrees instead of radians
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' anglecheck(c(0,1),radians=FALSE)
#'

anglecheck <- function(x,radians=TRUE){

if(x[1] == 0 & x[2] == 0) return(0)

theta <- atan(abs(x[2])/abs(x[1]))

theta2 <- (sign(x[1])==1)*(sign(x[2])==1)*theta +
		(sign(x[1])==-1)*(sign(x[2])==1)*(pi-theta) +
		(sign(x[1])==-1)*(sign(x[2])==-1)*(pi+theta) +
		(sign(x[1])==1)*(sign(x[2])==-1)*(2*pi-theta) +
		(x[1] == 0) * sign(x[2] > 0) * pi/2 +
		(x[1] == 0) * sign(x[2] < 0) * 3*pi/2 + 
		(x[1] > 0) * sign(x[2] == 0) * 0 +
		(x[1] < 0) * sign(x[2] == 0) * pi

if(radians == FALSE) {
return(theta2/2/pi * 360)
}
else return(theta2)

}
