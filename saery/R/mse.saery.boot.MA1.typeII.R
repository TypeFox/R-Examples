mse.saery.boot.MA1.typeII <-
function(X, D, md, beta, sigma2edi, sigmau1, sigmau2, theta, B){
  
  sigmau1.gorro.ast <- sigmau2.gorro.ast <- theta.gorro.ast <- vector()
  mudi.ast <- beta.gorro.ast <- mudi.gorro.ast <- u1di.gorro.ast <- u1di.gorro.ast2 <- u1d.gorro.ast <- u1d.gorro.ast2 <- u2di.gorro.ast <- list()
  u1 <- u1d <- u1di <- u2d <- u2di <- beta.gorro1.ast <- beta2.gorro.ast <- list()
  
  b <- BadTot2 <- 0
  M <- sum(md)
  mse.boot.MA1.typeII <- difmud <- 0
  excepcion <- vector()
  while(b<B){
    b <- b+1
    edi <- rnorm(M,0,sqrt(sigma2edi))
    u1di <- 0
    u2di <- 0
    for(d in 1:D){
      u1 <- rnorm(1,0,sqrt(sigmau1))
      u1d[[d]] <- rep(u1,md[d])
      epdi <- rnorm(md[d]+1,0,sqrt(sigmau2))
      for(i in 1:md[d])
        u2di[i] <- epdi[i+1] - theta*epdi[i]
      u2d[[d]] <- u2di
    }  
    u1di <- unlist(u1d)    
    u2di <- unlist(u2d)
    ydi.ast <- as.vector(X%*%beta + u1di + u2di + edi)                       
    mudi.ast[[b]] <- as.vector(X%*%beta + u1di + u2di)  
    sigma.0 <- sigmau1    
    
    fitboot <- try(REML.saery.MA1(X, ydi.ast, D, md, sigma2edi, sigma.0=sigma.0), TRUE)    # 
    
    if(class(fitboot)=="try-error"){
      excepcion <- c(excepcion, b)
      write.table(data.frame(class(fitboot),D,b,excepcion), file="WARNING.txt", append=TRUE, col.names=FALSE, row.names=FALSE)
      b <- b-1
    }
    else {
      cat("Bootstrap sample no. " , b,"\n")
      
      sigmau1.gorro.ast <- fitboot[[1]][1]
      sigmau2.gorro.ast <- fitboot[[1]][2]
      theta.gorro.ast <- fitboot[[1]][3]
      beta.gorro.ast[[b]] <- BETA.U.saery.MA1(X, ydi.ast, D, md, sigma2edi, sigmau1.gorro.ast, sigmau2.gorro.ast, theta.gorro.ast)   
      beta.gorro1.ast[[b]] <- beta.gorro.ast[[b]][[1]]
      u1d.gorro.ast[[b]] <- beta.gorro.ast[[b]][[2]]
      u1di.gorro.ast[[b]] <- list()
      for(d in 1:D)
        u1di.gorro.ast[[b]][[d]] <- rep(u1d.gorro.ast[[b]][d,1],md[d])
      u1di.gorro.ast[[b]] <- unlist(u1di.gorro.ast[[b]])
      u2di.gorro.ast[[b]] <- beta.gorro.ast[[b]][[3]]
      mudi.gorro.ast[[b]] <- as.vector(X %*%beta.gorro1.ast[[b]] + u1di.gorro.ast[[b]]+ u2di.gorro.ast[[b]])
      
      difmud <- difmud + (mudi.gorro.ast[[b]] - mudi.ast[[b]])^2
    }
  }
  mse.boot.MA1.typeII <- difmud/B
  
  return(mse.boot.MA1.typeII)
  
}
