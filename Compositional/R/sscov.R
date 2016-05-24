################################
#### Spatial sign covariance matrix
#### Tsagris Michail 2/2015  
#### References: A Durre, D Vogel, DE Tyler (2014)
#### The spatial sign covariance matrix with unknown location
#### Journal of Multivariate Analysis, 130: 107--117. 
#### http://arxiv.org/pdf/1307.5706v2.pdf 
#### mtsagris@yahoo.gr
################################

sscov <- function(x) {
  ## x contains the data
  x <- as.matrix(x)  ## makes sure x is a matrix
  n <- nrow(x)  ## sample size
  p <- ncol(x)
  me <- spat.med(x)  ## spatial median of x
  y <- x - rep( me, rep(n, p) )
  y <- y / sqrt ( rowSums(y^2) )  ## unit vectors
  crossprod( y ) / n  ## SSCM
} 