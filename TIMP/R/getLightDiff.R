"getLightDiff" <-
function(theta, x, kinscal, kmat, jvec,fixedkmat=FALSE,
         kinscalspecial = list(), kinscalspecialspec = list(),
         lightregimespec = list())
{
  depkIn <- lightregimespec$depkIn
  br <- lightregimespec$lightBreaks
  br <- append(br, tail(x,1) ) 
  inten <- lightregimespec$intensity   

  ## depkIn is 1 for those rates where there is intensity dependence 
  depkin <- as.logical(depkIn)
  
  dimk <- ncol(kmat)
    
  # get eigen for the case that the light is off
  theta1 <- theta
  theta1[which(depkin)] <- 0 
  
  K1 <- fillK(theta1, kinscal, kmat, fixedkmat, kinscalspecial,
             kinscalspecialspec,nocolsums=TRUE)
  eigenlijk1 <- eigen(K1, only.values = F)
  V1 <- eigenlijk1$vectors
  gamma1 <- solve(V1) 
  k1 <- -eigenlijk1$values
  cat("eigen dark", k1, "\n")
    
  # get eigen for the case that the light is on
  theta2 <- theta
  theta2[which(depkin)] <- theta[which(depkin)] * inten   
  K2 <- fillK(theta2, kinscal, kmat, fixedkmat, kinscalspecial,
             kinscalspecialspec, nocolsums=TRUE)
  eigenlijk2 <- eigen(K2, only.values = F)
  V2 <- eigenlijk2$vectors
  gamma2 <- solve(V2)
  k2 <- -eigenlijk2$values
  cat("eigen light", k2, "\n")
  
  A <- matrix(nrow = dimk, ncol = dimk)
  
  ## start with light off
  break2 <- 0
  
  for(i in 1:(length(br))) {
    
    break1 <- break2 + 1
    break2 <- which(x >= br[i])[1] 

    if(is.na(break2)) break2 <- length(x)
    if(break1 > break2)  break
    if(break2 < length(x)) break2 <- break2 - 1
    
    cat("times are from", break1, "to", break2,"\n")
    
    xtemp <- x[break1:break2]
    if(i %% 2 != 0) {
      k <- k1
      V <- V1
      gamma <- gamma1
    }
    else {
      k <- k2
      V <- V2
      gamma <- gamma2
    }

    if(i > 1) {
      # take the last concentrations
      finalcon <- c.tempA[nrow(c.tempA),] 
    } else finalcon <- jvec
    
    finalcon <- finalcon/max(finalcon)
    cat("finalcon", finalcon, "\n")
    gamma <- gamma %*% finalcon
    for(j in 1:dimk) 
      A[j, ] <- V[ ,j] * gamma[j]

    c.temp <- calcC(k, xtemp)
    c.temp <- c.temp %*% A
    
    if(i == 1)
      c.tempA <- c.temp
    else 
      c.tempA <- rbind(c.tempA, c.temp)
      
  }
  c.tempA
}
