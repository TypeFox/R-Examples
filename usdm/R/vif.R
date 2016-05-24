# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2013
# Version 1.1
# Licence GPL v3

.vif <- function(y) {
  z<-rep(NA,ncol(y))
  names(z) <- colnames(y)
  for (i in 1:ncol(y)) {
    z[i] <-  1/(1-summary(lm(y[,i]~.,data=y[-i]))$r.squared)
  }
  return(z)
}

.vif2 <- function(y,w) {
  z<-rep(NA,length(w))
  names(z) <- colnames(y)[w]
  for (i in 1:length(w)) {
    z[i] <-  1/(1-summary(lm(as.formula(paste(colnames(y)[w[i]],"~.",sep='')),data=y))$r.squared)
  }
  return(z)
}

.maxCor <- function(k){
  n <- nrow(k)
  for (i in 1:n) k[i:n,i] <- NA
  w <- which.max(k)
  c(rownames(k)[((w%/%nrow(k))+1)],colnames(k)[w%%nrow(k)])
}



.minCor <- function(k){
  k <- abs(k)
  rr<-c();cc<-c();co<-c()
  for (c in 1:(ncol(k)-1)) {
    for (r in (c+1):nrow(k)){
      rr<-c(rr,rownames(k)[r]);cc<-c(cc,colnames(k)[c])
      co <- c(co,k[r,c])
    }
  }
  w <- which.min(co)
  c(rr[w],cc[w])
}
#-----------------
if (!isGeneric("vif")) {
  setGeneric("vif", function(x, ...)
    standardGeneric("vif"))
}  
setMethod('vif', signature(x='RasterStackBrick'),
          function(x, maxobservations=5000) {
            if (nlayers(x) == 1) stop("The Raster object should have at least two layers")
            x <- as.data.frame(x)
            x <- na.omit(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            v <- .vif(x)
            data.frame(Variables=names(v),VIF=as.vector(v))
          }
)

setMethod('vif', signature(x='data.frame'),
          function(x, maxobservations=5000) {
            if (ncol(x) == 1) stop("At least two variables are needed to quantify vif")
            x <- na.omit(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            v <- .vif(x)
            data.frame(Variables=names(v),VIF=as.vector(v))
          }
)

setMethod('vif', signature(x='matrix'),
          function(x, maxobservations=5000) {
            if (ncol(x) == 1) stop("At least two variables are needed to quantify vif")
            x <- na.omit(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            v <- .vif(x)
            data.frame(Variables=names(v),VIF=as.vector(v))
          }
)

if (!isGeneric("vifcor")) {
  setGeneric("vifcor", function(x, th= 0.9, ...)
    standardGeneric("vifcor"))
}  
setMethod('vifcor', signature(x='RasterStackBrick'),
          function(x, th=0.9, maxobservations=5000) {
            if (nlayers(x) == 1) stop("The Raster object should have at least two layers")
            LOOP <- TRUE
            x <- as.data.frame(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            x <- na.omit(x)
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            while (LOOP) {
              xcor <- abs(cor(x))
              mx <- .maxCor(xcor)
              if (xcor[mx[1],mx[2]] >= th) {
                w1 <- which(colnames(xcor) == mx[1])
                w2 <- which(rownames(xcor) == mx[2])
                v <- .vif2(x,c(w1,w2))
                ex <- mx[which.max(v[mx])]
                exc <- c(exc,ex)
                x <- x[,-which(colnames(x) == ex)]
              } else LOOP <- FALSE
            }
            if (length(exc) > 0) n@excluded <- exc
            v <- .vif(x)
            n@corMatrix <- cor(x)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)

setMethod('vifcor', signature(x='data.frame'),
          function(x, th=0.9, maxobservations=5000) {
            if (ncol(x) == 1) stop("The Raster object should have at least two layers")
            LOOP <- TRUE
            x <- na.omit(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            while (LOOP) {
              xcor <- abs(cor(x))
              mx <- .maxCor(xcor)
              if (xcor[mx[1],mx[2]] >= th) {
                w1 <- which(colnames(xcor) == mx[1])
                w2 <- which(rownames(xcor) == mx[2])
                v <- .vif2(x,c(w1,w2))
                ex <- mx[which.max(v[mx])]
                exc <- c(exc,ex)
                x <- x[,-which(colnames(x) == ex)]
              } else LOOP <- FALSE
            }
            if (length(exc) > 0) n@excluded <- exc
            v <- .vif(x)
            n@corMatrix <- cor(x)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)

setMethod('vifcor', signature(x='matrix'),
          function(x, th=0.9, maxobservations=5000) {
            if (ncol(x) == 1) stop("The Raster object should have at least two layers")
            LOOP <- TRUE
            x <- as.data.frame(x)
            x <- na.omit(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            while (LOOP) {
              xcor <- abs(cor(x))
              mx <- .maxCor(xcor)
              if (xcor[mx[1],mx[2]] >= th) {
                w1 <- which(colnames(xcor) == mx[1])
                w2 <- which(rownames(xcor) == mx[2])
                v <- .vif2(x,c(w1,w2))
                ex <- mx[which.max(v[mx])]
                exc <- c(exc,ex)
                x <- x[,-which(colnames(x) == ex)]
              } else LOOP <- FALSE
            }
            if (length(exc) > 0) n@excluded <- exc
            v <- .vif(x)
            n@corMatrix <- cor(x)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)

if (!isGeneric("vifstep")) {
  setGeneric("vifstep", function(x, th= 10, ...)
    standardGeneric("vifstep"))
}

setMethod('vifstep', signature(x='RasterStackBrick'),
          function(x, th=10, maxobservations=5000) {
            if (nlayers(x) == 1) stop("The Raster object should have at least two layers")
            LOOP <- TRUE
            x <- as.data.frame(x)
            x <- na.omit(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            while (LOOP) {
              v <- .vif(x)
              if (v[which.max(v)] >= th) {
                ex <- names(v[which.max(v)])
                exc <- c(exc,ex)
                x <- x[,-which(colnames(x) == ex)]
              } else LOOP=FALSE
            }
            if (length(exc) > 0) n@excluded <- exc
            v <- .vif(x)
            n@corMatrix <- cor(x)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)

setMethod('vifstep', signature(x='data.frame'),
          function(x, th=10, maxobservations=5000) {
            if (ncol(x) == 1) stop("The Raster object should have at least two layers")
            LOOP <- TRUE
            x <- na.omit(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            while (LOOP) {
              v <- .vif(x)
              if (v[which.max(v)] >= th) {
                ex <- names(v[which.max(v)])
                exc <- c(exc,ex)
                x <- x[,-which(colnames(x) == ex)]
              } else LOOP=FALSE
            }
            if (length(exc) > 0) n@excluded <- exc
            v <- .vif(x)
            n@corMatrix <- cor(x)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)

setMethod('vifstep', signature(x='matrix'),
          function(x, th=10, maxobservations=5000) {
            if (ncol(x) == 1) stop("The Raster object should have at least two layers")
            LOOP <- TRUE
            x <- as.data.frame(x)
            x <- na.omit(x)
            if(nrow(x) > maxobservations) x <- x[sample(1:nrow(x),maxobservations),]
            n <- new("VIF")
            n@variables <- colnames(x)
            exc <- c()
            while (LOOP) {
              v <- .vif(x)
              if (v[which.max(v)] >= th) {
                ex <- names(v[which.max(v)])
                exc <- c(exc,ex)
                x <- x[,-which(colnames(x) == ex)]
              } else LOOP=FALSE
            }
            if (length(exc) > 0) n@excluded <- exc
            v <- .vif(x)
            n@corMatrix <- cor(x)
            n@results <- data.frame(Variables=names(v),VIF=as.vector(v))
            n
          }
)

