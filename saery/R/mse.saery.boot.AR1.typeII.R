mse.saery.boot.AR1.typeII <-
function(X, D, md, beta, sigma2edi, sigmau1, sigmau2, rho, B){
  
  sigmau1.gorro.ast <- sigmau2.gorro.ast <- rho.gorro.ast <-vector()
  u1<-u1d<-u1di<-mudi.ast.1 <-u2d<-u2di <- beta.gorro.ast<- beta.gorro1.ast<- mudi.gorro.ast.1 <- u1di.gorro.ast <- u1di.gorro.ast2<- u1d.gorro.ast <- u1d.gorro.ast2<-u2di.gorro.ast<-list()
  
  Badg2 <- 0
  b <- BadTot2 <- 0
  M <- sum(md)
  mse.boot.AR1.typeII <- difmud <- 0
  
  excepcion <- vector()
  while(b<B){
    b <- b+1
    edi <- rnorm(M,0,sqrt(sigma2edi))
    u1di <- 0
    u2di <- 0
    for(d in 1:D){
      u1 <- rnorm(1,0,sqrt(sigmau1))
      u1d[[d]] <- rep(u1,md[d])
      epdi <- rnorm(md[d],0,sqrt(sigmau2))
      u2di <- (1/(1-rho^2)^0.5)*epdi[1]
      for(i in 2:md[d])
        u2di[i] <- rho*u2di[i-1] + epdi[i]
      u2d[[d]] <- u2di
    }  
    u1di <- unlist(u1d)    
    u2di <- unlist(u2d)
    ydi.ast1 <- as.vector(X%*%beta + u1di + u2di + edi)                       
    mudi.ast.1[[b]] <- as.vector(X%*%beta + u1di + u2di)  
    sigma.0 <- sigmau1
    
    fitboot <- try(REML.saery.AR1(X, ydi.ast1, D, md, sigma2edi, sigma1.0=sigma.0,sigma2.0=sigma.0), TRUE)    # 
    
    if(class(fitboot)=="try-error"){
      excepcion <- c(excepcion, b)
      write.table(data.frame(class(fitboot),D,b), file="WARNING.txt", append=TRUE, col.names=FALSE, row.names=FALSE)
      b <- b-1
    }
    else {
      cat("Bootstrap sample no. " , b,"\n")
      
      sigmau1.gorro.ast <- fitboot[[1]][1]
      sigmau2.gorro.ast <- fitboot[[1]][2]
      rho.gorro.ast <- fitboot[[1]][3]
      beta.gorro.ast[[b]] <- BETA.U.saery.AR1(X, ydi.ast1, D, md, sigma2edi, sigmau1.gorro.ast, sigmau2.gorro.ast, rho.gorro.ast)
      beta.gorro1.ast[[b]] <- beta.gorro.ast[[b]][[1]]
      u1d.gorro.ast[[b]] <- beta.gorro.ast[[b]][[2]]
      u1di.gorro.ast[[b]] <- list()
      for(d in 1:D)
        u1di.gorro.ast[[b]][[d]] <- rep(u1d.gorro.ast[[b]][d,1],md[d])
      u1di.gorro.ast[[b]] <- unlist(u1di.gorro.ast[[b]])
      u2di.gorro.ast[[b]] <- beta.gorro.ast[[b]][[3]]
      mudi.gorro.ast.1[[b]] <- as.vector(X%*%beta.gorro1.ast[[b]] + u1di.gorro.ast[[b]]+ u2di.gorro.ast[[b]])
      
      difmud <- difmud+(mudi.gorro.ast.1[[b]]-mudi.ast.1[[b]])^2 
    }
  }
  mse.boot.AR1.typeII <- difmud/B
  
  return(mse.boot.AR1.typeII)
  
}
