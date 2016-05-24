#############################################################
#
#	sample.linearcomb function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: September, 1, 2008
#	Version: 0.2
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

#library(geometry)
#sample.linearcomb <- function(x, size) {
##  posch <- unique(as.vector(convhulln(data.matrix(x),"Qt")))
##  x <- data.matrix(x[posch,])
#  x <- data.matrix(x)
#  nn <- nrow(x)
#  res <- matrix(0, nrow=size, ncol=ncol(x))
#  for (i in 1:size) {
#    delta <- runif(nn, 0, 1)
#    delta <- delta/sum(delta)
#    res[i,] <- c(delta%*%x)
#  }
#  return(res)
#}

sample.linearcomb <- function(x, size, weights=rep(1, nrow(x)), subsize) {
  x <- data.matrix(x)
  nn <- nrow(x)
  res <- matrix(0, nrow=size, ncol=ncol(x))
  for (i in 1:size) {
    res[i,] <- apply(x[sample(1:nrow(x), size=subsize, replace=TRUE, prob=weights),], 2, mean)
  }
  return(res)
}
