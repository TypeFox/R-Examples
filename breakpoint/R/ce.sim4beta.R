ce.sim4beta <-
function(N, data, h, L0, L, M, Melite, eps, a){
  
  if (N == 0){
    seql <- c(0, L)
    mBic.full <- mBIC(seql, data, 0, L, h)    
    return(list(locis = c(1, (L + 1)), mBIC = mBic.full))
    rm(mBic.full, seql)
    
  } else {
    
  ########################Parameter initialization######################################################  
  new.para <- array(1, dim = c(2, N))
  ######################################################################################################  
#  modbic <- c()
  k <- 0
  repeat
  {
    k <- k + 1
    ch <- array(0, dim =c (M, (N + 2)))       
    ch[, 1] <- c(1)                    
    ch[, (N + 2)] <- c(L + 1)    
    ch[, (2 : (N + 1))] <- apply(new.para, 2, betarand, L0, L, M)   
    ch <- t(apply(ch, 1, sort))                       
    mod.bic <- apply(ch, 1, mBIC, data, N, L, h)
    ch <- cbind(ch, mod.bic)                         
    ch <- ch[order(ch[, (N + 3)], decreasing = TRUE), ]  
    melitesmpl <- ch[1 : Melite, ]                     
#    modbic[k] <- melitesmpl[1, (N + 3)] 
    
    newpar.n <- array(0, dim = c(2, N))
    newpar.n[1, ] <- apply(as.matrix(melitesmpl[, (2 :(N + 1))]), 2, mean)
    newpar.n[2, ] <- apply(as.matrix(melitesmpl[, (2 :(N + 1))]), 2, var)   
    
    newpar.new <- array(0, dim = c(2, N))
    newpar.new[1, ] <- apply(newpar.n, 2, fun.alpha, L0, L)
    newpar.new[2, ] <- apply(newpar.n, 2, fun.beta, L0, L)
    new.para <- a * newpar.new + (1 - a) * new.para
    
    mad <- apply(as.matrix(melitesmpl[, (2 : (N + 1))]), 2, mad)
    
    if(max(mad) <= eps){break}
  }
#  return(list(loci = ch[1, (1 : (N + 2))], mBIC = modbic[k]))
  return(list(loci = ch[1, (1 : (N + 2))], mBIC = melitesmpl[1, (N + 3)][[1]]))
  }
}
