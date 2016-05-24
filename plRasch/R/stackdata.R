### The function "plStackData" reads in the orginal examinee x items
### matrix, and generate data for pseudolikelihood estimation by using
### "the stacking data trick"
###
### it is a general function that can be used for multiple latent
### variables
###
### Arguments:
###   data : the orginal examinee response matrix (nexaminee x nitem)
###   item.mtx : adjacency matrix between items and latent traits
###   trait.mtx: adjacency matrix among latent traits
###
### Value:
###   A matrix of the stacked data is returned. The columns are
###   'resp': the response
###   'item': the item that is condictioned
###   'phi..': corresponding to the parameters in LMA
###   'cid': the strata id
###
### Notes:
###   for simplification, we use 'raw' scores for the repsonpse
### data, score matrix and data maxtrix are the same( 'data' in
### arguments)
###
### Author: Zhushan "Mandy" Li
### Last Updated: July 29, 2006

plStackData <- function(data, item.mtx, trait.mtx){

  ## make sure the arguments are matrices
  if(! is.matrix(data)) data <- as.matrix(data)
  if(! is.matrix(item.mtx)) item.mtx <- as.matrix(item.mtx)
  if(! is.matrix(trait.mtx)) trait.mtx <- as.matrix(trait.mtx)
  
  nexaminee <- nrow(data)  # number of examinees
  nitem <- ncol(data)  # number of items
  ntrait <- ncol(item.mtx)  # number of latent trait

  ## lower trianglular of the latent traits adjacency matrix
  trait.low <- trait.mtx & lower.tri(trait.mtx,diag=TRUE)
  ## total number of parameters for latent variables
  nphi <- sum(trait.low)

  which.low <- which(trait.low, arr.ind=TRUE)
  ## the names of the parameters
  phi.names <- paste("phi",apply(which.low, 1, paste, collapse=''), sep='')

  ## construct matricies used for calculating phi
  phi.array <- array(0, dim=c(ntrait, ntrait, nphi))
  phi.array[cbind(which.low, 1:nphi)] <- 1
  phi.array[cbind(which.low[,2:1], 1:nphi)] <- 1
  
  jjmtx <- array(dim=c(nitem, nitem, nphi))
  for (i in 1:nphi){
    jjmtx[,,i] <- item.mtx %*% phi.array[,,i] %*% t(item.mtx)
    diag(jjmtx[,,i]) <- 0
  }
  
  ### stacking the data
  ## caluculate phi
  phi <- matrix(NA, nrow=nitem*nexaminee, ncol=nphi)
  colnames(phi) <- phi.names
  for(i in 1:nphi){
    phi[,i] <- as.vector(data %*% jjmtx[,,i])
  }

  ## the response 
  resp <- as.vector(data)
  ## item
  item <- rep(1:nitem, each=nexaminee)
  ## strata id. the rows in the new data genereated a sigle row in old
  ## data will have the same id 
  cid <- rep(1:nexaminee, nitem); 

  newdata <- cbind(resp, item, phi, cid)
  return(newdata);
}

