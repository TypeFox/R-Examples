ce.simNormalNB.AIC.Init <-
function(N, init.locs, data, h, L0, L, M, Melite, eps, a, b, r, var.init){
  
########################Parameter initialization######################################################  
    new_para <- rbind(init.locs, rep(var.init, N))    
######################################################################################################  
    
#    aic <- c()
    k <- 0
    
    repeat
    {
      k <- k + 1
      ch <- array(0, dim = c(M, N + 2))       
      ch[, 1] <- c(1)                    
      ch[, ( N + 2)] <- c(L + 1)    
      ch[, (2:(N + 1))] <- apply(new.para, 2, normrand, L0, L, M)       
      ch <- t(apply(ch, 1, sort))                              
      loglike <- apply(ch, 1, llhoodnb, data, r, h)
      aic.vals <- apply(as.data.frame(loglike), 1, AICnb, N)       
      ch <- cbind(ch, aic.vals, loglike)
      ch <- ch[order(ch[, (N + 3)], decreasing = FALSE), ]
      melitesmpl <- ch[1:Melite, ]                     
#      aic[k] <- melitesmpl[1, (N + 3)] 
      
      newpar.n <- array(0, dim = c(2, N))
      newpar.n[1, ] <-apply(as.matrix(melitesmpl[, (2:(N + 1))]), 2, mean)
      newpar.n[2, ] <-apply(as.matrix(melitesmpl[, (2:(N + 1))]), 2, sd)   
      
      new.para[1, ] <- a * newpar.n[1, ] + (1 - a) * new.para[1, ]
      new.para[2, ] <- b * newpar.n[2, ] + (1 - b) * new.para[2, ]
      
      mad <- apply(as.matrix(melitesmpl[, (2:(N + 1))]), 2, mad)
      
      if(max(mad) <= eps){break}
    }
#    return(list(loci = ch[1, (1:(N + 2))], BIC = bic[k]))
    return(list(loci = ch[1, (1:(N + 2))], AIC = melitesmpl[1, (N + 3)][[1]], LogLike = melitesmpl[1, (N + 4)][[1]]))
}
