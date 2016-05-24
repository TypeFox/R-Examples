ce.simNormalZINB <-
function(N, data, h, L0, L, M, Melite, eps, a, b, r){
  
  if (N == 0){
    loglik.full <- loglikzinb(1, (L+1), data, r, h)[[1]]
    BIC.full <- BICzinb(loglik.full, 0, L)
    return(list(loci = c(1, (L + 1)), BIC = BIC.full))
    rm(BIC.full, loglik.full)
    
  } else {
    
########################Parameter initialization######################################################  
    new.para <- rbind(rep(L0 + (L - L0)/2, N), rep(sqrt(L - L0)^2/12), N)    
######################################################################################################  
    
    bic <- c()
    k <- 0
    
    repeat
    {
      k <- k + 1
      ch <- array(0, dim = c(M, N + 2))       
      ch[, 1] <- c(1)                    
      ch[, ( N + 2)] <- c(L + 1)    
      ch[, (2:(N + 1))] <- apply(new.para, 2, normrand, L0, L, M)       
      ch <- t(apply(ch, 1, sort))                              
      loglike <- apply(ch, 1, llhoodzinb, data, r, h)
      bic.vals <- apply(as.data.frame(loglike), 1, BICzinb, N, L)       
      ch <- cbind(ch, bic.vals)
      ch <- ch[order(ch[, (N + 3)], decreasing = FALSE), ]
      melitesmpl <- ch[1:Melite, ]                     
      bic[k] <- melitesmpl[1, (N + 3)] 
      
      newpar.n <- array(0, dim = c(2, N))
      newpar.n[1, ] <-apply(as.matrix(melitesmpl[, (2:(N + 1))]), 2, mean)
      newpar.n[2, ] <-apply(as.matrix(melitesmpl[, (2:(N + 1))]), 2, sd)   
      
      new.para[1, ] <- a * newpar.n[1, ] + (1 - a) * new.para[1, ]
      new.para[2, ] <- b * newpar.n[2, ] + (1 - b) * new.para[2, ]
      
      mad <- apply(as.matrix(melitesmpl[, (2:(N + 1))]), 2, mad)
      
      if(max(mad) <= eps){break}
    }
    return(list(loci = ch[1, (1:(N + 2))], BIC = bic[k]))
  }
}
