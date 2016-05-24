basePlot <-
function (K) {
## BASEPLOT Plot a contour of the 2D Gaussian distribution with covariance matrix K.
## FORMAT
## DESC Creates the basic plot as an ellipse with major and minor radii as the
## square root of the two eigenvalues respectively.
## ARG K : the covariance matrix.
##
## COPYRIGHT : Neil D. Lawrence, 2008; Alfredo Kalaitzis, 2010

  eigVecs = eigen(K)
  U = eigVecs$vectors[ , 2:1] ## Reverse order of eigenvectors(columns).
  V = eigVecs$values
  V[V<0] = as.complex(V[V<0])
  r = Re(sqrt(V))
  theta = seq(0, 2*pi, length=200)
  xy = cbind(r[1]*sin(theta), r[2]*cos(theta))%*%U
  plot(0, type = "n", xlim=c(min(xy[,1]), max(xy[,1])),
    ylim=c(min(xy[,2]), max(xy[,2])), xlab='', ylab='') 
  cont = lines(xy[, 1], xy[, 2], col='blue') ## 'lines' only applies on existing plots.

  zeroAxes()
}
