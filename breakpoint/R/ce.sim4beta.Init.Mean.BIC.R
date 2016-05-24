ce.sim4beta.Init.Mean.BIC <-
function(N, init.locs, data, h, L0, L, M, Melite, eps, a, var.init){
 
    v <- var(data[, 1])
    
#     if (N == 0){
#       LL.full <- apply(ch, 1, loglikMeanVarNormal, data, h)
#       BIC.val <- -2*LL.full + 2*(N + 1)* log(L)
#       #AIC.val <- -2*LL.full + 4* (N + 1)
#       return(list(locis = c(1, (L + 1)), BIC.Val = BIC.val, LogLike = LL.full))
#       rm(LL.full, BIC.val)
#       
#     } else {
#       
      ########################Parameter initialization######################################################  
      new.para <- array(0, dim = c(2, N))
      new.para[1, ] <- betaIntEst(init.locs[1:N], var.init, L0, L)$alpha.init
      new.para[2, ] <- betaIntEst(init.locs[1:N], var.init, L0, L)$beta.init
      ######################################################################################################  
#      llVal <- c()
#      bic <- c()
     # aic <- c()
      k <- 0
      repeat
      {
        k <- k + 1
        ch <- array(0, dim =c (M, (N + 2)))       
        ch[, 1] <- c(1)                    
        ch[, (N + 2)] <- c(L + 1)    
        ch[, (2 : (N + 1))] <- apply(new.para, 2, betarand, L0, L, M)   
        ch <- t(apply(ch, 1, sort))                       
        
        
#         LL.full <- apply(ch, 1, loglikMeanNormal, data, h)
#         BIC.val <- -2*LL.full + (N + 2)* log(L)
       # AIC.val <- -2*LL.full + 4* (N + 1)
        LL.full <- apply(ch, 1, llhood.MeanNormal, data, v, h)
        BIC.val <- apply(as.data.frame(LL.full), 1, BIC.MeanNormal, N, L) 
        
      #  mod.bic <- apply(ch, 1, mBIC, data, N, L, h)
      #  ch <- cbind(ch, mod.bic)
        ch <- cbind(ch, LL.full, BIC.val)
        ch <- ch[order(ch[, (N + 4)], decreasing = FALSE), ]  
        melitesmpl <- ch[1 : Melite, ]                     
        #modbic[k] <- melitesmpl[1, (N + 3)] 
#        llVal[k] <- melitesmpl[1, (N + 3)] 
#        bic[k] <- melitesmpl[1, (N + 4)] 
        
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
#      return(list(loci = ch[1, (1 : (N + 2))], BIC.Val = bic[k], LogLike = llVal[k]))
      return(list(loci = ch[1, (1 : (N + 2))], BIC.Val = melitesmpl[1, (N + 4)][[1]], LogLike = melitesmpl[1, (N + 3)][[1]]))
  }
