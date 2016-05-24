################################
#### Rotation matrix 
#### Tsagris Michail 10/2013 
#### mtsagris@yahoo.gr
#### References: G. J. A. Amaral, I. L. Dryden & Andrew T. A. Wood (2007, JASA)
#### Pivotal Bootstrap Methods for k-Sample Problems in Directional Statistics and Shape Analysis. 
################################

rotation <- function(a, b) {
  ## a and b are two unit vectors
  ## Calculates the rotation matrix
  ## to move a to b along the geodesic path
  ## on the unit sphere which connects a to b
  p <- length(a)
  ca <- a - b * crossprod(a, b)
  ca <- ca /sqrt( sum(ca^2) )
  A <- tcrossprod(b, ca) 
  A <- A - t(A)
  theta <- acos( sum(a * b) )
  diag(p) + sin(theta) * A + (cos(theta) - 1) * ( tcrossprod(b) + 
  tcrossprod(ca) )
}