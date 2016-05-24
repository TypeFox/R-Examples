#############################################################
#
#	localdepth2Dsimplicial function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: June, 15, 2008
#	Version: 0.5
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth2Dsimplicial <- function(x, y, tau, use) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  nx <- nrow(x)
  ny <- nrow(y)
  result <- list()

  ## funzione per calcolare la lunghezza di un vettore e l'area di un triangolo#
##  norm <- function(x) sqrt(t(x)%*%x)
##  norma <- function(x) sqrt(diag(x%*%t(x)))
##  area <- function(x) sqrt(sum(x)*(sum(x)/2-x[1])*(sum(x)/2-x[2])*(sum(x)/2-x[3])/2)
  ## numero di terne
  nt <- choose(nx, 3)
  depth <- rep(0,ny)
  localdepth <- rep(0,ny)
  if (use=='diameter') {
    res <- .C("ld2Ddiamsimp",
              x1 = as.double(x[,1]),
              x2 = as.double(x[,2]),
              y1 = as.double(y[,1]),
              y2 = as.double(y[,2]),
              tau =as.double(tau),
              nx = as.integer(nx),
              ny = as.integer(ny),
              depth = as.double(depth),
              localdepth = as.double(localdepth),
              DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")
  result$localdepth <- res$localdepth  
  result$depth <- res$depth
  result$num <- nt
  return(result)
  } else {
    res <- .C("ld2Dareasimp",
              x1 = as.double(x[,1]),
              x2 = as.double(x[,2]),
              y1 = as.double(y[,1]),
              y2 = as.double(y[,2]),
              tau =as.double(tau),
              nx = as.integer(nx),
              ny = as.integer(ny),
              depth = as.double(depth),
              localdepth = as.double(localdepth),
              DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")    
#   aree <- posizione <- rep(0, nt)
#   for (i1 in 1:(nx-2)) {
#   for (i2 in (i1+1):(nx-1)) {
#      for (i3 in (i2+1):nx) {
#       ind1 <- ind1+1
#        x1 <- x[i1,]
#        x2 <- x[i2,]
#        x3 <- x[i3,]
#        n1 <- norm(x1-x2)
#        n2 <- norm(x1-x3)
#        n3 <- norm(x2-x3)
#        aree[ind1] <- area(c(n1,n2,n3))
#      }
#    }
#  }
#  x <- t(x)
#  for (ind2 in 1:ny) {
#    y1 <- y[ind2,]
#    u <- (x-y1)/norma(x-y1)
#    alpha <- atan2(u[2,],u[1,])%%(2*pi)
#    ind1 <- 0
#    for (i1 in 1:(nx-2)) {
#      for (i2 in (i1+1):(nx-1)) {
#        for (i3 in (i2+1):nx) {
#          ind1 <- ind1 + 1
#          a <- sort(alpha[c(i1,i2,i3)])
#          ## Ho tolto a[1]==0 ||
#          if ((a[3]-a[1])%%(2*pi) >= pi & (a[1]-a[2])%%(2*pi) >= pi & (a[2]-a[3])%%(2*pi) >= pi) {
#            depth[ind2] <- depth[ind2]+1
#            if (use=='areas' && aree[ind1] <= tau)
#              localdepth[ind2] <- localdepth[ind2]+1      
#            if (use=='spherical' && max(norma(x[,c(i1,i2,i3)]-y1))<= tau)
#              localdepth[ind2] <- localdepth[ind2]+1
#          }
#        }
#      }
#    }
#  }
  result$localdepth <- res$localdepth  
  result$depth <- res$depth
#  result$areas <- aree
  result$num <- nt
  return(result)
 }
}

