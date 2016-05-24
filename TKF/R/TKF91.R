#########################################################################
# File Name: ../R/TKF91.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Sun Sep  7 14:23:14 2014
#########################################################################

TKF91Pair <- function(seq1, seq2, mu=NULL, distance=NULL,
                  ## mu: by default is 0.001 from median of mu values from Fungi dataset.
                  expectedLength=362, 
                  substModel, substModelBF){
  if(!all(seq1 %in% AACharacterSet) || !all(seq2 %in% AACharacterSet)){
    stop("This implementation currently only supports 20 AA characters ", 
         paste(AACharacterSet, collapse=" "))
  }
  seq1Int <- AAToInt(seq1)
  seq2Int <- AAToInt(seq2)
  ## for the C matrix index
  seq1Int <- seq1Int - 1L
  seq2Int <- seq2Int - 1L

  expectedLength <- as.numeric(expectedLength)
  
  if(is.null(mu) && is.null(distance)){ 
    ## Do the 2D optimisation
    ans <- .Call("TKF91LikelihoodFunction2DMainNM", seq1Int, seq2Int,
                 expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, expectedLength, substModel, substModelBF){
                          ansTemp <- .Call("TKF91LikelihoodFunctionWrapper", seq1Int, seq2Int, x[1], x[2], expectedLength, substModel, substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, c(ans["PAM"], ans["Mu"]),
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 expectedLength=expectedLength, substModel=substModel,
                 substModelBF=substModelBF)
    return(c(ans, "PAMVariance"=solve(ansHessian)[1,1],
             "MuVariance"=solve(ansHessian)[2,2],
             "coVariance"=solve(ansHessian)[1,2]))
  }else if(!is.null(mu) && is.null(distance)){
    ## Do the 1D distance optimisation
    ans <- .Call("TKF91LikelihoodFunction1DMain", seq1Int, seq2Int, mu,
                 expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, mu, expectedLength, substModel, substModelBF){
                          ansTemp <- .Call("TKF91LikelihoodFunctionWrapper", seq1Int, seq2Int, x, mu, expectedLength, substModel, substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, ans["PAM"], 
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 mu=mu, expectedLength=expectedLength,
                 substModel=substModel, 
                 substModelBF=substModelBF)
    return(c(ans, "PAMVariance"=solve(ansHessian)[1,1]))
  }else if(!is.null(mu) && !is.null(distance)){
    ## Just calculate the likelihood, given mu and distance
    ans <- .Call("TKF91LikelihoodFunctionWrapper", seq1Int, seq2Int, 
                 distance, mu, expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, expectedLength, substModel, substModelBF){
                          ansTemp <- .Call("TKF91LikelihoodFunctionWrapper", seq1Int, seq2Int, x[1], x[2], expectedLength, substModel, substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, c(ans["PAM"], ans["Mu"]),
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 expectedLength=expectedLength, substModel=substModel,
                 substModelBF=substModelBF)
    return(c(ans, "PAMVariance"=solve(ansHessian)[1,1],
             "MuVariance"=solve(ansHessian)[2,2],
             "coVariance"=solve(ansHessian)[1,2]))
  }else{
    stop("You cannot estimate mu alone!")
  }
}


TKF91 <- function(fasta, mu=NULL, expectedLength=362, 
                  substModel, substModelBF){
  seqnames <- names(fasta)
  nSeqs <- length(fasta)
  distanceMatrix <- matrix(0, ncol=nSeqs, nrow=nSeqs,
                           dimnames=list(seqnames, seqnames))
  varianceMatrix <- distanceMatrix
  if(is.null(mu)){
    muMatrix <- distanceMatrix
  }
  negLoglikelihoodMatrix <- distanceMatrix  
  for(i in 1:(nSeqs-1L)){
    for(j in (i+1L):nSeqs){
      message(seqnames[i], " vs ", seqnames[j])
      ans <- TKF91Pair(fasta[[i]], fasta[[j]], 
                       mu=mu, expectedLength=expectedLength,
                       substModel=substModel, substModelBF=substModelBF)
      distanceMatrix[i,j] <- distanceMatrix[j,i] <- ans["PAM"]
      varianceMatrix[i,j] <- varianceMatrix[j,i] <- ans["PAMVariance"]
      negLoglikelihoodMatrix[i,j] <- negLoglikelihoodMatrix[j,i] <- 
        ans["negLogLikelihood"]
      if(is.null(mu)){
        muMatrix[i,j] <- muMatrix[j,i] <- ans["Mu"]
      }
    }
  }
  if(is.null(mu)){
    return(list(distanceMatrix=distanceMatrix,
                varianceMatrix=varianceMatrix,
                muMatrix=muMatrix,
                negLoglikelihoodMatrix=negLoglikelihoodMatrix
                ))
  }else{
    return(list(distanceMatrix=distanceMatrix,
                varianceMatrix=varianceMatrix,
                negLoglikelihoodMatrix=negLoglikelihoodMatrix
                ))
  }
}



