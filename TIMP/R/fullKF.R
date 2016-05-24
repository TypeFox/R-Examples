"fullKF" <-
  function(theta, kinscal, kmat, jvec,fixedkmat=FALSE,
           kinscalspecial = list(), kinscalspecialspec = list(),
           nocolsums = FALSE)
  {
    #fullKF
    
    dimk <- nrow(kmat)
    A <- matrix(nrow = dimk, ncol = dimk)
    K <- fillK(theta, kinscal, kmat, fixedkmat, kinscalspecial,
               kinscalspecialspec, nocolsums)
    eigenlijk <- eigen(K, only.values = F)
    V <- eigenlijk$vectors
    gamma <- solve(V) %*% jvec
    for(j in 1:dimk) 
      A[j, ] <- V[ ,j] * gamma[j]
    list(A=A, values = eigenlijk$values)
  }
