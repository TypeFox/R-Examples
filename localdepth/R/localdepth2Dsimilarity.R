#############################################################
#
#	localdepth2Dsimplicialsimilarity function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: October, 02, 2008
#	Version: 0.2-1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth2Dsimplicialsimilarity <- function(x, y, tau, use, weight=NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  nx <- nrow(x)
  ny <- nrow(y)
  result <- list()

  ## numero di terne
  nt <- choose(nx, 3)
  depth <- rep(0, ny*nt)
  localdepth <- rep(0, ny*nt)
  if (use=='diameter') {
    res <- .C("ld2Ddiamsimpsim",
              x1 = as.double(x[,1]),
              x2 = as.double(x[,2]),
              y1 = as.double(y[,1]),
              y2 = as.double(y[,2]),
              tau =as.double(tau),
              nx = as.integer(nx),
              ny = as.integer(ny),
              nt = as.integer(nt),
              depth = as.double(depth),
              localdepth = as.double(localdepth),
              DUP = TRUE, NAOK = FALSE, PACKAGE = "localdepth")
  } else {
    res <- .C("ld2Dareasimpsim",
              x1 = as.double(x[,1]),
              x2 = as.double(x[,2]),
              y1 = as.double(y[,1]),
              y2 = as.double(y[,2]),
              tau =as.double(tau),
              nx = as.integer(nx),
              ny = as.integer(ny),
              nt = as.integer(nt),              
              depth = as.double(depth),
              localdepth = as.double(localdepth),
              DUP = TRUE, NAOK = FALSE, PACKAGE = "localdepth")
  }
  res$depth <- matrix(res$depth, nrow=ny, ncol=nt, byrow=FALSE)
  res$localdepth <- matrix(res$localdepth, nrow=ny, ncol=nt, byrow=FALSE)
  
#  tetracorica <- function(x, y) {
#    a <- x%*%y
#    b <- x%*%(1-y)
#    cc <- (1-x)%*%y
#    d <- (1-x)%*%(1-y)
#    a <- c(a, b, cc, d)
#    return(a)
#  }
#  result$jaccard <- result$czekanowski <- result$sokal <- matrix(0, ny, ny)
#  for (i in 1:ny) {
#    for (j in i:ny) {
#      a <- tetracorica(res$localdepth[i,], res$localdepth[j,])
#      result$jaccard[i,j] <- result$jaccard[j,i] <- a[1]/(a[1]+a[2]+a[3])
#      result$czekanowski[i,j] <- result$czekanowski[j,i] <- 2*a[1]/(2*a[1]+a[2]+a[3])
#      result$sokal[i,j] <- result$sokal[j,i] <- a[1]+a[4]/sum(a)
#    }
#  }
##  result$ld <- apply(res$localdepth, 1, sum)/nt
##  result$d <- apply(res$depth, 1, sum)/nt
##  result$ld <- res$localdepth
##  result$d <- res$depth
  if (!is.null(weight)) {
    result <- rep(0, nt)
    if (use=='diameter') {
      vol <- .C("twoDdiam", x = as.double(x[,1]), y = as.double(x[,2]),
               nx = as.integer(nx), result = as.double(result),
               DUP = TRUE, NAOK = FALSE, PACKAGE = "localdepth")$result
    } else {
      vol <- .C("twoDarea", x = as.double(x[,1]), y = as.double(x[,2]),
               nx = as.integer(nx), result = as.double(result),
               DUP = TRUE, NAOK = FALSE, PACKAGE = "localdepth")$result
    }
    vol <- weight(vol)
    res$depth <- res$depth%*%(t(res$depth)*vol)
    res$localdepth <- res$localdepth%*%(t(res$localdepth)*vol)
  } else {
    res$depth <- res$depth%*%t(res$depth)
    res$localdepth <- res$localdepth%*%t(res$localdepth)
  }

  result$localdepth <- res$localdepth
  result$depth <- res$depth
  result$call <- match.call()
  result$num <- nt
  class(result) <- 'localdepth.similarity'
  return(result)
}
