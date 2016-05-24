#' Computes Multiplier Bootstrap realizations of the supremum of a Gaussian 
#' field on a contour.
#'
#' @param x x-Coordinates of the grid on which the data is observed.
#' @param y y-Coordinates of the grid on which the data is observed.
#' @param R An array of dimension c(length(x),length(y),n) containing the 
#'          realizations of the field.
#' @param f A matrix of dimension c(length(x),length(y)). The contour on which 
#'          tail probabilities are to be computed is defind as {f=c}. 
#' @param level The level of the contour.
#' @param N The number of Bootstrap realizations to produce. Default is 1000.

#' @return A vector of length N containing the Bootstrap realizations of the 
#'         supremum. 
MBContour = function(x, y, R, f, level, N=1000){

  cont = contourLines(list(x=x,y=y,z=f),levels=level,nlevels=1)
  if(length(cont)==0) return(rep(-Inf,N))
  cont_x = c()
  cont_y = c()
  for(l in 1:length(cont)) {cont_x = c(cont_x,cont[[l]]$x); cont_y = c(cont_y,cont[[l]]$y)}
  cont = cbind(cont_x,cont_y)
  
  interp_max = function(G){
    G = matrix(G,ncol=length(y))
    max(fields::interp.surface(list(x=x,y=y,z=G),cont))
  }
  
  n = dim(R)[3]
  R = matrix(R,ncol=n)
  g = matrix(rnorm(n*N),n,N)
  apply(abs(R %*% g ),2,interp_max) / sqrt(n-2)
}