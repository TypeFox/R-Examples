#' This function plots the weights interacting between estimated effects for each pixel. 
#'
#' 
#'
#' @name sim.weightsplot
#' @aliases sim.weightsplot
#' @title Plot Function for Weights (Simulated Data)
#' @usage sim.weightsplot(weights, nei, nx, ny, thresh = 0.05, ...) 
#' @param weights matrix, containing MCMC-output the of posterior estimates of weights.
#' @param nei matrix, locations of weights in precision matrix.
#' @param nx scalar, number of pixels in x-direction.
#' @param ny scalar, number of pixels in y-direction.
#' @param thresh scalar, defining the threshold to which the median of the weights smaller than 
#'        this threshold should be plotted.
#' @param \dots graphical parameters for \code{image} can also be passed on as arguments to 
#'       this function. 
#' @author Max Hughes
#' @note This function is solely for MCMC-outputs on simulated data.


sim.weightsplot <- function(weights, nei, nx, ny, thresh=0.05, ...) { 

  require(spatstat)

  #if(class(model)!="adpatsmoFMRI")
  #  stop("Inserted model needs to be of class spatsmoFMRI.")

  if(is.matrix(weights)!=TRUE)
    stop("weights needs to be a matrix!")
  if(dim(nei)[2]!=dim(weights)[1])
    stop("The row dimension of weights needs to be the same size as the
          column dimension of nei!")
  
  nrow <- nx
  ncol <- ny
  w <- apply(weights,1,median)
  
  if(ncol(nei)!=((ncol-1)*nrow + (nrow-1)*ncol))
    stop("Wrong matrix dimension!")
  
  getcoor <- function(image, valbool) {
    if(is.logical(valbool)) 
      ind <- which(valbool)
    else
      ind <- which(image$v==valbool)
    nx <- ncol(image$v)
    ny <- nrow(image$v)
    xind <- ceiling(ind/ny)
    X <- image$xcol[xind]
    Y <- image$yrow[ind - (xind-1)*ny]
    data.frame(X,Y)
  }
  
  nrowim <- 2*nrow - 1
  ncolim <- 2*ncol - 1

  Z <- matrix(nrow=nrowim, ncol=ncolim)
  Zim <- im(Z, xcol=1:ncolim, yrow=1:nrowim)

  inds <- rep(seq(1,ncolim,by=2), times=nrow) + rep((0:(nrow-1))*(2*ncolim), each=ncol)

  pixelinds <- rep(0, nrowim*ncolim)

  pixelinds[inds] <- 1:(nrow*ncol)

  Zim$v <- matrix(ncol=ncolim, nrow=nrowim, data=pixelinds, byrow=TRUE)

  coorweights <- sapply(as.data.frame(nei), function(x) colMeans(getcoor(Zim, Zim$v %in% x)))

  coorweights <- coorweights[,w<thresh]
  
  Zim[list(x=coorweights[1,], y=coorweights[2,])] <- 1.5

  coorzero <- getcoor(Zim, Zim$v==0)
  Zim[list(x=coorzero$X, y=coorzero$Y)] <- NA

  coornonweights <- getcoor(Zim, Zim$v!=1.5)
  Zim[list(x=coornonweights$X, y=coornonweights$Y)] <- 1

  coorweights <- getcoor(Zim, Zim$v==1.5)
  Zim[list(x=coorweights$X, y=coorweights$Y)] <- 2

  Zim$v <- Zim$v[1:nrow(Zim$v),1:ncol(Zim$v)]   
  image(Zim$v, x=1:nrow(Zim$v), y=1:ncol(Zim$v), axes=FALSE, xlab="", ylab="", ...)
}

