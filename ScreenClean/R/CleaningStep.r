
#####################################################
#
#  This file contains the following functions
#
#  VectorizeBase(): express an integer on some base as a vector
#
#  PMLE(): fit the penalized MLE 
#
#  CleaningStep(): do the cleaning step of the GS/UPS procedures.
#  
####################################################

VectorizeBase <- function(i, base, length){
  #############################################################################
  #
  # express the number i on the base as a vector.
  #
  # Args: 
  #   i: the nonnegative number to be converted.
  #   base: the base of the convertion.
  #   length: the length of the output vector
  #
  # Returns: 
  #   vectorizeBase: the converted vector
  #
  #############################################################################
  
  if ((i < 0) || (floor(i) != i)) {
  stop("i must be a non-negative integer.")
  }
  vector <- rep(0, length)
  for( j in length:1) {
    vector[j] <- i - base * floor(i/base)
    i <- (i - vector[j]) / base;
  }
  return(vector)
}

PMLE <- function(gram, y, lambda, uu) {
  ###########################################################################
  #
  # fits the penalized MLE in the cleaning step of the GS/UPS procedures.
  # It requires quadprog package.
  #
  # Args:
  #   gram: the sub gram matrix of the small scale problem.
  #   y: the sub-vector of y.tilde 
  #   lambda, uu: the tuning parameters of GS, they have the intuitive 
  #              interpretation of the desired sparse level and
  #              the minimal signal strength to be detected. 
  #
  ###########################################################################
  
  n <- length(y);
  b <- matrix(0,nrow=n, ncol=1)
  l <- 0  
  # run through all possible signs 
  for( k in 0:(3^n-1)){
    # express k on base 3 as a vector of length n
    idx <- VectorizeBase(k, 3, n)
    # convert to a col of -1, 0, 1 instead of 0,1,2
    idx <- matrix(idx,ncol=1) - 1; 
    cluster <- (idx != 0)
    card <- sum(cluster)
    if(card == 0){
      lt <- 0
      bt <- matrix(0, nrow=n, ncol=1)
    }else{
      signs <- idx[cluster]
      yc <- y[cluster];
      gramcluster <- as.matrix(as.matrix(gram)[cluster,cluster])
      uucluster <- uu * rep(1, card) 
      # expand the signs to a matrix
      amat <- diag(card)
      diag(amat) <- signs;
      ###The following function is from quadprog package. Check its help document for more details. 
      minpoint <- solve.QP(gramcluster, yc, amat, uucluster)      
      # note constrOptim is not reliable for one
      # solve.QP is faster than nloptr
      bt <- minpoint$solution
      lt <- minpoint$value + lambda^2 * card/2  # add the constant part
      if(lt<l){
        b <- matrix(0, nrow=n, ncol=1)
        b[cluster] <- bt;
        l <- lt
      }    
    }      
  }
return(b)  
}
  

CleaningStep <- function(survivor, y.tilde, gram, lambda, uu){  
  ############################################################################# 
  #
  # runs the cleaning step of the GS/UPS procedure.
  #
  # Args: 
  #   survivor: the result of the screening step, a logical vector.
  #   y.tilde: y.tilde = X'%*%y, where X and y are the predictor matrix and the reponse.
  #   gram: the thresholded sparse gram matrix.
  #   lambda, uu: the tuning parameters of GS, they have the intuitive 
  #            interpretation of the desired sparse level and the minimal
  #            signal strength to be detected. 
  #
  # Returns:
  #   beta.gs: the estimated regression coefficient of the graphlet screening
  #
  #############################################################################
  gram<-as.matrix(gram)
  p <- dim(gram)[2]
  ##survivor <- as.logical(survivor)
  n.survivor <- sum(survivor)
  
  # extract for survivors
  yt <- y.tilde[survivor] # shorter
  omega <- as.matrix(gram[survivor, survivor])
  beta <- matrix(0, nrow=n.survivor, ncol=1) # col vector
  remain <- matrix(TRUE, nrow=n.survivor, ncol=1) # col vector 
  
  # find the maximum connected subgraphs 
  while(sum(remain) > 0) {
    i <- min(which(remain > .5))
    cluster <- matrix(FALSE, nrow=n.survivor, ncol=1)
    cluster[i] <- TRUE
    new.cluster <- matrix((omega[, i] != 0), ncol=1)
    while(!identical(cluster, new.cluster)) {
      cluster <- new.cluster
      new.cluster <- matrix((rowSums(abs(omega[, cluster])) != 0), ncol=1)
    }    
    # run penalized MLE for this cluster
    if(sum(cluster) > 20){
        print(paste("The cluster length is ",sum(cluster)))
      stop("cluster too long. The program has stopped.")
    }
    beta[cluster] <- PMLE(omega[cluster, cluster], yt[cluster], lambda, uu)
    remain[cluster] <- FALSE
  }  
  # expand the original dimension by adding zeros
  beta.gs <- matrix(0,nrow=p,ncol=1)
  beta.gs[survivor] <- beta

return(beta.gs)
}
  
