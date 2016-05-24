# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2013
# Version 1.2
# Licence GPL v3



.neighborRowCol <- function(x,cell,w) {
  nc <- ncol(x)
  nr <- nrow(x)
  row <- trunc((cell-1)/nc) + 1
  col <- cell - ((row-1) * nc)
  
  mr <- matrix(ncol=w,nrow=w) 
  tr <- trunc(w/2)
  
  for (i in 1:w) {
    mr[i,] <-tr 
    tr <- tr - 1
  }
  i <- trunc(w/2)+1
  mr[i,i] <- NA
  row <- row + as(mr,"vector")
  col <- col + as(-t(mr),"vector")
  row[row < 1 | row > nr] <- NA
  col[col < 1 | col > nc] <- NA
  rc <- cbind(row, col)
  rc
}

.filter<-function(r,d1=0,d2) {
  c<- d2%/%r
  x<- y<- seq(-c,c,1)
  eg<- expand.grid(x,y)
  eg[1]<- -eg[1]
  eg[3]<- sqrt(eg[1]^2+eg[2]^2)
  ndim<- c*2+1
  m<-matrix(eg[,3],ncol=ndim,nrow=ndim)*r
  mw<- which(m > d1 & m <= d2)
  m[mw] = 1 
  m[-mw] = NA
  mw <- trunc(length(m)/2)+1
  m[mw]<- 1
  return (m)
}


.localK<-function(x,d1=0,d2,cell) {
  filter<-.filter(r=res(x)[1],d1=d1,d2=d2)
  if (length(filter) < 3) stop("d2 is less than resolution!")
  xv <- x[1:ncell(x)]
  w <- which(!is.na(xv))
  xv <- xv[w]
  ras <- FALSE
  if (missing(cell)) {
    cell <-c(1:ncell(x))[w]
    ras <- TRUE
    xx <- raster(x)
  }
  nf <- ncol(filter)
  n <- length(xv)
  n1 <- n-1
  n2 <- n-2
  
  x <- as(x,"matrix")
  out <- rep(NA,length(cell))
  for (c in 1:length(cell)) {
    xi <-  t(x)[cell[c]]
    if (!is.na(xi)) {
      xn <- .neighborRowCol(x,cell[c],nf)
      xn[,1] <- xn[,1] * as(filter,"vector")
      xn <- na.omit(xn)
      xn <- unlist(lapply(1:nrow(xn),function(r) {x[xn[r,1],xn[r,2]]}))
      xn <- xn[!is.na(xn)]
      zj <- xv - xi
      Di <- sum(abs(zj))/n1
      Wi <- length(xn)
      Ki <- sum(abs(xi-xn))
      Fi <- sum(zj*zj)/n1
      VKi <- Wi*(n1-Wi) * (Fi-Di*Di) / n2
      out[c] <- (Ki - Wi*Di) / sqrt(VKi)
    } else out[c] <- NA
  }
  if (ras) {
    xx[cell] <- out
    return(xx)
  } else return(out)
}




.localGeary <- function(x,d1=0,d2,cell) {
  filter<-.filter(r=res(x)[1],d1=d1,d2=d2)
  if (length(filter) < 3) stop("d2 is less than resolution!")
  sigmaX <- cellStats(x,sum)
  i <- trunc(length(filter)/2)+1
  ras <- FALSE
  if (missing(cell)) {
    w <- which(!is.na(x[1:ncell(x)]))
    n <- length(w)
    cell <-c(1:ncell(x))[w]
    rm(w)
    ras <- TRUE
    xx <- raster(x)
  } else n <- ncell(x) - cellStats(x,"countNA")
  nf <- ncol(filter)
  n1 <- n-1
  s2 <- cellStats(x,var)
  s2 <- (s2 * n1)/n
  x <- as(x,"matrix")
  out <- rep(NA,length(cell))
  for (c in 1:length(cell)) {
    xi <-  t(x)[cell[c]]
    if (!is.na(xi)) {
      xn <- .neighborRowCol(x,cell[c],nf)
      xn[,1] <- xn[,1] * as(filter,"vector")
      xn <- na.omit(xn)
      xn <- unlist(lapply(1:nrow(xn),function(r) {x[xn[r,1],xn[r,2]]}))
      xn <- xn[!is.na(xn)]
      Eij <- sum((xi - xn)^2)
      out[c] <- Eij / s2
      
    } else out[c] <- NA
  }
  
  if (ras) {
    xx[cell] <- out
    return(xx)
  } else return(out)
}

.localG <- function(x,d1=0,d2,cell) {
  filter<-.filter(r=res(x)[1],d1=d1,d2=d2)
  if (length(filter) < 3) stop("d2 is less than resolution!")
  sigmaX <- cellStats(x,sum)
  sigmaX2 <- cellStats(x^2,sum)
  i <- trunc(length(filter)/2)+1
  ras <- FALSE
  if (missing(cell)) {
    w <- which(!is.na(x[1:ncell(x)]))
    n <- length(w)
    cell <-c(1:ncell(x))[w]
    rm(w)
    ras <- TRUE
    xx <- raster(x)
  } else n <- ncell(x) - cellStats(x,"countNA")
  nf <- ncol(filter)
  n1 <- n - 1
  n2 <- n - 2
  s2 <- cellStats(x,var)
  s2 <- (s2 * n1)/n
  x <- as(x,"matrix")
  out <- rep(NA,length(cell))
  for (c in 1:length(cell)) {
    xi <-  t(x)[cell[c]]
    if (!is.na(xi)) {
      xn <- .neighborRowCol(x,cell[c],nf)
      xn[,1] <- xn[,1] * as(filter,"vector")
      xn <- na.omit(xn)
      xn <- unlist(lapply(1:nrow(xn),function(r) {x[xn[r,1],xn[r,2]]}))
      xn <- xn[!is.na(xn)]
      s.xi <- sum(xn)
      xbar.i <- (sigmaX - xi)/n1
      Wi <- length(xn)
      si <- sqrt(((sigmaX2 - xi^2) / n1) - xbar.i^2)
      G <- (s.xi - Wi * xbar.i) / (si * sqrt(((n1*Wi) - Wi^2) / n2))
      out[c] <- G
      
    } else out[c] <- NA
  }
  
  if (ras) {
    xx[cell] <- out
    return(xx)
  } else return(out)
}


.localMoran <- function(x,d1=0,d2,cell) {
  filter<-.filter(r=res(x)[1],d1=d1,d2=d2)
  if (length(filter) < 3) stop("d2 is less than resolution!")
  xv <- x[1:ncell(x)]
  w <- which(!is.na(xv))
  xv <- xv[w]
  i <- trunc(length(filter)/2)+1
  ras <- FALSE
  if (missing(cell)) {
    cell <-c(1:ncell(x))[w]
    ras <- TRUE
    xx <- raster(x)
    } 
  nf <- ncol(filter)
  n <- length(xv)
  n1 <- n-1
  n2 <- n-2
  sigmaX <- sum(xv)
  sigmaX2 <- sum(xv * xv)
  
  s2 <- var(xv)
  s2 <- (s2 * n1)/n
  s4 <- sum((xv - mean(xv))^4)/n
  rm(xv,w)
  b2 <- s4 / (s2 ^ 2)
  x <- as(x,"matrix")
  out <- rep(NA,length(cell))
  for (c in 1:length(cell)) {
    xi <-  t(x)[cell[c]]
    if (!is.na(xi)) {
      xn <- .neighborRowCol(x,cell[c],nf)
      xn[,1] <- xn[,1] * as(filter,"vector")
      xn <- na.omit(xn)
      xn <- unlist(lapply(1:nrow(xn),function(r) {x[xn[r,1],xn[r,2]]}))
      xn <- xn[!is.na(xn)]
      s.x <- (sigmaX-xi)
      xbar= s.x / n1
      s <- ((sigmaX2 - xi * xi) - (s.x * s.x / n1)) / n1
      z <- xn- xbar
      lz <- sum(z)
      Ii <- ((xi-xbar) / s) * lz
      Wi <- length(xn)
      EI <- -Wi/n1
      VarI <- (Wi * (n - b2) / n1) + (Wi^2 * ((2*b2 - n) / (n1*n2))) - (Wi / n1^2)
      out[c] <- (Ii - EI) / sqrt(VarI)
    } else out[c] <- NA
  }
  
  if (ras) {
    xx[cell] <- out
    return(xx)
  } else return(out)
}

.localG2 <- function(x,d1=0,d2,cell) {
  filter<-.filter(r=res(x)[1],d1=d1,d2=d2)
  if (length(filter) < 3) stop("d2 is less than resolution!")
  sigmaX <- cellStats(x,sum)
  sigmaX2 <- cellStats(x^2,sum)
  i <- trunc(length(filter)/2)+1
  ras <- FALSE
  if (missing(cell)) {
    w <- which(!is.na(x[1:ncell(x)]))
    n <- length(w)
    cell <-c(1:ncell(x))[w]
    rm(w)
    ras <- TRUE
    xx <- raster(x)
  } else n <- ncell(x) - cellStats(x,"countNA")
  nf <- ncol(filter)
  n1 <- n - 1
  n2 <- n - 2
  xbar <- sigmaX / n
  s2 <- (sigmaX2 / n) - (xbar ^ 2)
  si <- sqrt(s2)
  x <- as(x,"matrix")
  out <- rep(NA,length(cell))
  for (c in 1:length(cell)) {
    xi <-  t(x)[cell[c]]
    if (!is.na(xi)) {
      xn <- .neighborRowCol(x,cell[c],nf)
      xn[,1] <- xn[,1] * as(filter,"vector")
      xn <- na.omit(xn)
      xn <- unlist(lapply(1:nrow(xn),function(r) {x[xn[r,1],xn[r,2]]}))
      xn <- xn[!is.na(xn)]
      s.xi <- sum(xn) + xi
      
      Wi <- length(xn)+1
      
      G <- (s.xi - Wi * xbar) / (si * sqrt(((n*Wi) - Wi^2) / n1))
      out[c] <- G
      
    } else out[c] <- NA
  }
  
  if (ras) {
    xx[cell] <- out
    return(xx)
  } else return(out)
}


if (!isGeneric("lisa")) {
  setGeneric("lisa", function(x, y, d1=0, d2, cell, statistic="I")
    standardGeneric("lisa"))
}  

setMethod('lisa', signature(x='Raster',y='missing'), 
          function(x, y, d1=0, d2, cell, statistic="I") {
            if (!statistic %in% c("K1","I","G","G*","C","k1","i","c","g","g*")) stop("statistic should be one of K1, I, G, G*, and C")
            if (nlayers(x) == 1) {
              if (missing(cell)) {
                if (statistic == "K1" | statistic == "k1") return(.localK(x,d1=d1,d2=d2))
                if (statistic == "I" | statistic == "i") return(.localMoran(x,d1=d1,d2=d2))
                if (statistic == "G" | statistic == "g") return(.localG(x,d1=d1,d2=d2))
                if (statistic == "G*" | statistic == "g*") return(.localG2(x,d1=d1,d2=d2))
                if (statistic == "C" | statistic == "c") return(.localGeary(x,d1=d1,d2=d2))
              } else {
                if (statistic == "K1" | statistic == "k1") return(.localK(x,d1=d1,d2=d2,cell=cell))
                if (statistic == "I" | statistic == "i") return(.localMoran(x,d1=d1,d2=d2,cell=cell))
                if (statistic == "G" | statistic == "g") return(.localG(x,d1=d1,d2=d2,cell=cell))
                if (statistic == "G*" | statistic == "g*") return(.localG2(x,d1=d1,d2=d2,cell=cell))
                if (statistic == "C" | statistic == "c") return(.localGeary(x,d1=d1,d2=d2,cell=cell))
              }
            } else {
              if (missing(cell)) {
                out <- raster(x[[1]])
                if (statistic == "K1" | statistic == "k1") {
                  for (i in 1:nlayers(x)) out <- addLayer(out,.localK(x[[i]],d1=d1,d2=d2))
                  names(out) <- paste(statistic," _statistic_for_",names(x),sep="")
                  out <- brick(out)
                }
                if (statistic == "I" | statistic == "i") { 
                  for (i in 1:nlayers(x)) out <- addLayer(out,.localMoran(x[[i]],d1=d1,d2=d2))
                  names(out) <- paste(statistic," _statistic_for_",names(x),sep="")
                  out <- brick(out)
                }
                if (statistic == "G" | statistic == "g") {
                  for (i in 1:nlayers(x)) out <- addLayer(out,.localG(x[[i]],d1=d1,d2=d2))
                  names(out) <- paste(statistic," _statistic_for_",names(x),sep="")
                  out <- brick(out)
                }
                if (statistic == "G*" | statistic == "g*") {
                  for (i in 1:nlayers(x)) out <- addLayer(out,.localG2(x[[i]],d1=d1,d2=d2))
                  names(out) <- paste(statistic," _statistic_for_",names(x),sep="")
                  out <- brick(out)
                }
                if (statistic == "C" | statistic == "c") {
                  for (i in 1:nlayers(x)) out <- addLayer(out,.localGeary(x[[i]],d1=d1,d2=d2))
                  names(out) <- paste(statistic," _statistic_for_",names(x),sep="")
                  out <- brick(out)
                }
              } else {
                out <- matrix(nrow=length(cell),ncol=nlayers(x))
                colnames(out) <- names(x)
                rownames(out) <- cell
                if (statistic == "K1" | statistic == "k1") for (i in 1:nlayers(x)) out[,i] <- .localK(x[[i]],d1=d1,d2=d2,cell=cell)
                if (statistic == "I" | statistic == "i") for (i in 1:nlayers(x)) out[,i] <- .localMoran(x[[i]],d1=d1,d2=d2,cell=cell)
                if (statistic == "G" | statistic == "g") for (i in 1:nlayers(x)) out[,i] <- .localG(x[[i]],d1=d1,d2=d2,cell=cell)
                if (statistic == "G*" | statistic == "g*") for (i in 1:nlayers(x)) out[,i] <- .localG2(x[[i]],d1=d1,d2=d2,cell=cell)
                if (statistic == "C" | statistic == "c") for (i in 1:nlayers(x)) out[,i] <- .localGeary(x[[i]],d1=d1,d2=d2,cell=cell)
              }
              out
            }
            
          }
)

setMethod('lisa', signature(x='Raster',y='SpatialPoints'), 
          function(x, y, d1=0, d2, cell, statistic="I") {
            if (!statistic %in% c("K1","I","G","G*","C")) stop("statistic should be one of K1, I, G, G*, and C")
            
            xy <- coordinates(y)
            cell <- cellFromXY(x,xy)
            
            if (nlayers(x) == 1) {
                if (statistic == "K1" | statistic == "k1") return(.localK(x,d1=d1,d2=d2,cell=cell))
                if (statistic == "I" | statistic == "i") return(.localMoran(x,d1=d1,d2=d2,cell=cell))
                if (statistic == "G" | statistic == "g") return(.localG(x,d1=d1,d2=d2,cell=cell))
                if (statistic == "G*" | statistic == "g*") return(.localG2(x,d1=d1,d2=d2,cell=cell))
                if (statistic == "C" | statistic == "c") return(.localGeary(x,d1=d1,d2=d2,cell=cell))
              } else {
                out <- matrix(nrow=length(cell),ncol=nlayers(x))
                colnames(out) <- names(x)
                rownames(out) <- cell
                if (statistic == "K1" | statistic == "k1") for (i in 1:nlayers(x)) out[,i] <- .localK(x[[i]],d1=d1,d2=d2,cell=cell)
                if (statistic == "I" | statistic == "i") for (i in 1:nlayers(x)) out[,i] <- .localMoran(x[[i]],d1=d1,d2=d2,cell=cell)
                if (statistic == "G" | statistic == "g") for (i in 1:nlayers(x)) out[,i] <- .localG(x[[i]],d1=d1,d2=d2,cell=cell)
                if (statistic == "G*" | statistic == "g*") for (i in 1:nlayers(x)) out[,i] <- .localG2(x[[i]],d1=d1,d2=d2,cell=cell)
                if (statistic == "C" | statistic == "c") for (i in 1:nlayers(x)) out[,i] <- .localGeary(x[[i]],d1=d1,d2=d2,cell=cell)
                out
              }
           }
)

setMethod('lisa', signature(x='Raster',y='SpatialPointsDataFrame'), 
          function(x, y, d1=0, d2, cell, statistic="I") {
            if (!statistic %in% c("K1","I","G","G*","C")) stop("statistic should be one of K1, I, G, G*, and C")
            
            xy <- coordinates(y)
            cell <- cellFromXY(x,xy)
            
            if (nlayers(x) == 1) {
              if (statistic == "K1" | statistic == "k1") return(.localK(x,d1=d1,d2=d2,cell=cell))
              if (statistic == "I" | statistic == "i") return(.localMoran(x,d1=d1,d2=d2,cell=cell))
              if (statistic == "G" | statistic == "g") return(.localG(x,d1=d1,d2=d2,cell=cell))
              if (statistic == "G*" | statistic == "g*") return(.localG2(x,d1=d1,d2=d2,cell=cell))
              if (statistic == "C" | statistic == "c") return(.localGeary(x,d1=d1,d2=d2,cell=cell))
            } else {
              out <- matrix(nrow=length(cell),ncol=nlayers(x))
              colnames(out) <- names(x)
              rownames(out) <- cell
              if (statistic == "K1" | statistic == "k1") for (i in 1:nlayers(x)) out[,i] <- .localK(x[[i]],d1=d1,d2=d2,cell=cell)
              if (statistic == "I" | statistic == "i") for (i in 1:nlayers(x)) out[,i] <- .localMoran(x[[i]],d1=d1,d2=d2,cell=cell)
              if (statistic == "G" | statistic == "g") for (i in 1:nlayers(x)) out[,i] <- .localG(x[[i]],d1=d1,d2=d2,cell=cell)
              if (statistic == "G*" | statistic == "g*") for (i in 1:nlayers(x)) out[,i] <- .localG2(x[[i]],d1=d1,d2=d2,cell=cell)
              if (statistic == "C" | statistic == "c") for (i in 1:nlayers(x)) out[,i] <- .localGeary(x[[i]],d1=d1,d2=d2,cell=cell)
              out
            }
            
          }
)