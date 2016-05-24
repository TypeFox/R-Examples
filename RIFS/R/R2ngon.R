# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function:
# R2ngon() - a regular polygonal protofractal set in R^2.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# n1, n2 - numbers of vertices & partition points for 
#          the edges of a regular polygon; 
# r, o - the radius & center of the circumscribed circle;
# cycle - logical; if cycle=FALSE, first & last points are not equal.
# Variables:
# phi - a vector of angular coordinates of all polygon vertices.
# Value:
# a matrix of points coordinates of a protofractal set in R^2.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
R2ngon <- function(n1=3, n2=1, r=1, o=c(0,0), cycle=FALSE) {
  phi <- seq(0, 2*pi, length=n1+1)
  x <- approx(o[1]+r*cos(phi), n=n1*n2+1)$y
  y <- approx(o[2]+r*sin(phi), n=n1*n2+1)$y
  if (cycle) return(cbind(x=x,y=y))
  else return(cbind(x=x[-1],y=y[-1]))
}