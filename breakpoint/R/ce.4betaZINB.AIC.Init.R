ce.4betaZINB.AIC.Init <-
function(N, init.locs, data, h, L0, L, M, Melite, eps, a, r, var.init){
  
  ######################### Parameter initialization #####################
  new.para <- array(0, dim = c(2, N))
  new.para[1, ] <- betaIntEst(init.locs[1:N], var.init, L0, L)$alpha.init
  new.para[2, ] <- betaIntEst(init.locs[1:N], var.init, L0, L)$beta.init
  ########################################################################
  
   k <- 0
  
  repeat
  {
    k <- k + 1
    ch <- array(0, dim = c(M, (N + 2)))       
    ch[, 1] <- c(1)                    
    ch[, (N + 2)] <- c(L + 1)    
    ch[, (2 : (N + 1))] <- apply(new.para, 2, betarand, L0, L, M)   
    ch <- t(apply(ch, 1, sort))                              
    loglike <- apply(ch, 1, llhoodzinb, data, r, h)
    aic.vals <- apply(as.data.frame(loglike), 1, AICzinb, N, L) 
    ch <- cbind(ch, aic.vals, loglike)
    ch <- ch[order(ch[, (N + 3)], decreasing = FALSE), ]
    melitesmpl <- ch[1:Melite, ]                     
#    aic[k] <- melitesmpl[1, (N + 3)] 

    newpar.n <- array(0, dim = c(2, N))
    newpar.n[1, ] <- apply(as.matrix(melitesmpl[, (2 : (N + 1))]), 2, mean)
    newpar.n[2, ] <- apply(as.matrix(melitesmpl[, (2 : (N + 1))]), 2, var)   
    
    newpar.new <- array(0, dim = c(2, N))
    newpar.new[1, ] <- apply(newpar.n, 2, fun.alpha, L0, L)
    newpar.new[2, ] <- apply(newpar.n, 2, fun.beta, L0, L)
    new.para <- a * newpar.new + (1 - a) * new.para
    
    mad <- apply(as.matrix(melitesmpl[, (2 : (N + 1))]), 2, mad)
    
    if(max(mad) <= eps){break}
  }
    # return(list(loci = ch[1, (1 : (N + 2))], AIC = aic[k]))
  return(list(loci = ch[1, (1 : (N + 2))], AIC = melitesmpl[1, (N + 3)][[1]], LogLike = melitesmpl[1, (N + 4)][[1]]))
  }
