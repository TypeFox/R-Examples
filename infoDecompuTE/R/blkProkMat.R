#pre- and post-multiply NTginvATN by block projection matrices for each stratum 
blkProkMat <- 
  function(z, T, N) {
    if (!is.matrix(z)) 
        return(z)
    
    nEffect = length(T)
    PNTginvATNP = T
    
    PNTginvATNP[[1]] = z %*% N %*% T[[1]] %*% invInfMat(C = z, N = N, T = T[[1]]) %*% T[[1]] %*% t(N) %*% t(z)
    
    newZ = (z %*% t(z)) - PNTginvATNP[[1]]
    
    if (nEffect != 1) {
        for (i in 2:nEffect) {
            PNTginvATNP[[i]] = newZ %*% N %*% t(T[[i]]) %*% invInfMat(C = newZ, N = N, T = T[[i]]) %*% T[[i]] %*% t(N) %*% t(newZ)
            newZ = (newZ %*% t(newZ)) - PNTginvATNP[[i]] 
        }
    }
    
    PNTginvATNP$Residual = newZ
    elementToRemove = numeric(0)
    for (i in 1:length(PNTginvATNP)) {
        if (all(PNTginvATNP[[i]] < 1e-06)) 
            elementToRemove = c(elementToRemove, i)
    }
    
    if (length(elementToRemove) > 0) 
        PNTginvATNP = PNTginvATNP[-elementToRemove]
    
    return(PNTginvATNP)
} 

