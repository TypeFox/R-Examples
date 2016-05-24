mse.saery.boot.MA1.typeI <-
function(X, D, md, beta, sigma2edi, sigmau1, sigmau2, theta, Fsig, B){
  
  sigmau1.gorro.ast <- sigmau2.gorro.ast <- theta.gorro.ast <- vector()
  u1g <- u1gdi <- u1 <- u1d <- mudi.ast <- u2d <- beta.gorro.ast <- mudi.gorro.ast <- u1di.gorro.ast <- u1di.gorro.ast2 <- u1d.gorro.ast <- u1d.gorro.ast2 <- u2di.gorro.ast <- mudi.gorro.ast <- beta.tilde.ast <- mudi.tilde.ast <- u1di.tilde.ast <- u1di.tilde.ast2 <- u1d.tilde.ast <- u1d.tilde.ast2 <- u2di.tilde.ast <- beta.gorro1.ast <- beta.tilde1.ast <- list()
  u1 <- u1d <- u2d <- msedi_ast <- list()
  
  g1g2g3 <- g123.MA1(X,D,md,sigma2edi,sigmau1, sigmau2, theta, Fsig)
  mse.boot.MA1.typeI <- g1g2g3[[1]] + g1g2g3[[2]]
  
  Badg2 <- 0
  b <- BadTot2 <- 0
  difmud2 <- 0
  M <- sum(md)
  excepcion <- vector()
  while(b<B){
    b <- b+1
    edi <- rnorm(M,0,sqrt(sigma2edi))
    u1di <- 0
    u2di <- 0
    for(d in 1:D){
      u1 <- rnorm(1,0,sqrt(sigmau1))
      for(i in 1:md[d])
        u1di[i] <- u1
      u1d[[d]] <- u1di
      epdi <- rnorm(md[d]+1,0,sqrt(sigmau2))
      for(i in 1:md[d])
        u2di[i] <- epdi[i+1]-theta*epdi[i]
      u2d[[d]] <- u2di
    }  
    u1di <- unlist(u1d)    
    u2di <- unlist(u2d)
    ydi.ast <- X%*%beta + u1di + u2di + edi                       
    mudi.ast[[b]] <- X%*%beta + u1di + u2di
    sigma.0 <- sigmau1
    
    fit2 <- try(REML.saery.MA1(X, ydi.ast, D, md, sigma2edi, sigma.0=sigma.0, MAXITER=50), TRUE)    # 
    
    if(class(fit2)=="try-error"){
      excepcion <- c(excepcion, b)
      write.table(data.frame(class(fit2),D,b), file="WARNING.txt", append=TRUE, col.names=FALSE)
      b <- b-1
      BadTot2 <- BadTot2 + 1
    }
    else {
      cat("Bootstrap Sample no. " , b,"\n")
      
      if(fit2[[3]]<50) {  
        sigmau1.gorro.ast <- fit2[[1]][1]
        sigmau2.gorro.ast <- fit2[[1]][2]
        theta.gorro.ast <- fit2[[1]][3]
        beta.gorro.ast[[b]] <- BETA.U.saery.MA1(X, ydi.ast, D, md, sigma2edi, sigmau1.gorro.ast, sigmau2.gorro.ast, theta.gorro.ast)
        beta.gorro1.ast[[b]] <- beta.gorro.ast[[b]][[1]]
        u1d.gorro.ast[[b]] <- beta.gorro.ast[[b]][[2]]
        u1di.gorro.ast[[b]] <- list()
        for(d in 1:D)
          u1di.gorro.ast[[b]][[d]] <- rep(u1d.gorro.ast[[b]][d,1],md[d])
        u1di.gorro.ast[[b]] <- unlist(u1di.gorro.ast[[b]])
        u2di.gorro.ast[[b]] <- beta.gorro.ast[[b]][[3]]
        mudi.gorro.ast[[b]] <- as.vector(X%*%beta.gorro1.ast[[b]] + u1di.gorro.ast[[b]] + u2di.gorro.ast[[b]])
        
        beta.tilde.ast[[b]] <- BETA.U.saery.MA1(X, ydi.ast, D, md, sigma2edi, sigmau1, sigmau2, theta)   
        beta.tilde1.ast[[b]] <- beta.tilde.ast[[b]][[1]]
        u1d.tilde.ast[[b]] <- beta.tilde.ast[[b]][[2]]
        u1di.tilde.ast[[b]] <- list()
        for(d in 1:D)
          u1di.tilde.ast[[b]][[d]] <- rep(u1d.tilde.ast[[b]][d,1],md[d])
        u1di.tilde.ast[[b]] <- unlist(u1di.tilde.ast[[b]])
        u2di.tilde.ast[[b]] <- beta.tilde.ast[[b]][[3]]
        mudi.tilde.ast[[b]] <- as.vector(X%*%beta.tilde1.ast[[b]] + u1di.tilde.ast[[b]] + u2di.tilde.ast[[b]])
        
        difmud2 <- difmud2 + (mudi.gorro.ast[[b]]-mudi.tilde.ast[[b]])^2
      }
      else {
        b <- b-1
      }
    }
  }
  mse.boot.MA1.typeI <- mse.boot.MA1.typeI + difmud2/B
  
  return(mse.boot.MA1.typeI)
  
}
