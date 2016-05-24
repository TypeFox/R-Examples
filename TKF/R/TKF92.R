TKF92Pair <- function(seq1, seq2, mu=NULL, r=NULL, distance=NULL,
                      ## mu: by default is 0.001 from median of mu values 
                      ## of Fungi dataset.
                      ## r: the probability in the geometric distribution 
                      ## of fragment length. By default, it should be 0.8480
                      ## from median of r values of Fungi dataset.
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
  
  if(is.null(mu) && is.null(distance) && is.null(r)){ 
    ## Do the 3D optimisation
    ans <- .Call("TKF92LikelihoodFunction3DMainNM", seq1Int, seq2Int,
                 expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, expectedLength, substModel, substModelBF){
                          ansTemp <- .Call("TKF92LikelihoodFunctionWrapper", seq1Int, seq2Int, x[1], x[2], x[3], expectedLength, substModel, substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, c(ans["PAM"], ans["Mu"], ans["r"]),
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 expectedLength=expectedLength, substModel=substModel,
                 substModelBF=substModelBF)
    return(c(ans, "PAMVariance"=solve(ansHessian)[1,1],
             "MuVariance"=solve(ansHessian)[2,2],
             "rVariance"=solve(ansHessian)[3,3],
             "coVariancePAMMu"=solve(ansHessian)[1,2],
             "coVariancePAMr"=solve(ansHessian)[1,3],
             "coVarianceMur"=solve(ansHessian)[2,3]
             )
           )
  }else if(!is.null(mu) && is.null(distance) && !is.null(r)){
    ## Do the 1D distance optimisation
    ans <- .Call("TKF92LikelihoodFunction1DMain", seq1Int, seq2Int, mu, r, 
                 expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, mu, r, expectedLength, substModel, substModelBF){
                          ansTemp <- .Call("TKF92LikelihoodFunctionWrapper", seq1Int, seq2Int, x, mu, r, expectedLength, substModel, substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, ans["PAM"], 
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 mu=mu, r=r, expectedLength=expectedLength,
                 substModel=substModel, 
                 substModelBF=substModelBF)
    return(c(ans, "PAMVariance"=solve(ansHessian)[1,1]))
  }else if(!is.null(mu) && !is.null(distance) && !is.null(r)){
    ## Just calculate the likelihood, given mu and distance
    ans <- .Call("TKF92LikelihoodFunctionWrapper", seq1Int, seq2Int, 
                 distance, mu, r, expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, expectedLength, substModel, substModelBF){
                          ansTemp <- .Call("TKF92LikelihoodFunctionWrapper", seq1Int, seq2Int, x[1], x[2], x[3], expectedLength, substModel, substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, c(ans["PAM"], ans["Mu"], ans["r"]),
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 expectedLength=expectedLength, substModel=substModel,
                 substModelBF=substModelBF)
    return(c(ans, "PAMVariance"=solve(ansHessian)[1,1],
             "MuVariance"=solve(ansHessian)[2,2],
             "rVariance"=solve(ansHessian)[3,3],
             "coVariancePAMMu"=solve(ansHessian)[1,2],
             "coVariancePAMr"=solve(ansHessian)[1,3],
             "coVarianceMur"=solve(ansHessian)[2,3]
             )
           )
  }else{
    stop("You cannot estimate mu or r alone!")
  }
}

TKF92 <- function(fasta, mu=NULL, r=NULL, expectedLength=362,
                  substModel, substModelBF){
  seqnames <- names(fasta)
  nSeqs <- length(fasta)
  distanceMatrix <- matrix(0, ncol=nSeqs, nrow=nSeqs,
                           dimnames=list(seqnames, seqnames))
  varianceMatrix <- distanceMatrix
  negLoglikelihoodMatrix <- distanceMatrix
  if(is.null(mu) && is.null(r)){
    muMatrix <- distanceMatrix
    rMatrix <- distanceMatrix
  }
  for(i in 1:(nSeqs-1L)){
    for(j in (i+1L):nSeqs){
      message(seqnames[i], " vs ", seqnames[j])
      ans <- TKF92Pair(fasta[[i]], fasta[[j]],
                       mu=mu, r=r, expectedLength=expectedLength,
                       substModel=substModel, substModelBF=substModelBF)
      distanceMatrix[i,j] <- distanceMatrix[j,i] <- ans["PAM"]
      varianceMatrix[i,j] <- varianceMatrix[j,i] <- ans["PAMVariance"]
      negLoglikelihoodMatrix[i,j] <- negLoglikelihoodMatrix[j,i] <-
        ans["negLogLikelihood"]
      if(is.null(mu) && is.null(r)){
        muMatrix[i,j] <- muMatrix[j,i] <- ans["Mu"]
        rMatrix[i,j] <- rMatrix[j,i] <- ans["r"]
      }
    }
  }
  if(is.null(mu) && is.null(r)){
    return(list(distanceMatrix=distanceMatrix,
                varianceMatrix=varianceMatrix,
                muMatrix=muMatrix,
                rMatrix=rMatrix,
                negLoglikelihoodMatrix=negLoglikelihoodMatrix))
  }else{
    return(list(distanceMatrix=distanceMatrix,
                varianceMatrix=varianceMatrix,
                negLoglikelihoodMatrix=negLoglikelihoodMatrix))
  }
}



