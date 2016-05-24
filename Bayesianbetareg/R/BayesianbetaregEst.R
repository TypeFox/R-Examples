BayesianbetaregEst <-
function (Y,X,Z,nsim,bpri,Bpri,gpri,Gpri,burn,jump,graph1,graph2){
  
  
 
  
  ### Cepeda - Metropolis - Hastings
  
  Y=as.matrix(Y)
  
  if(is.null(X)|is.null(Z)|is.null(Y)){
    stop("There is no data")
  }
  
  if(burn> 1 | burn < 0){
    stop("Burn must be a proportion between 0 and 1")
  }
  
  if(nsim <= 0){
    stop("the number of simulations must be greater than 0")
  }
  
  if(jump < 0|jump > nsim){
    stop("Jumper must be a positive number lesser than nsim")
  }
  
  
  
  ind1<-rep(0,nsim)
  ind2<-rep(0,nsim)
  betas.ind <- matrix(bpri,nrow=ncol(X))
  gammas.ind <- matrix(gpri,nrow=ncol(Z))
  
  beta.mcmc<-matrix(NA,nrow=nsim,ncol=ncol(X))
  gamma.mcmc<-matrix(NA,nrow=nsim,ncol=ncol(Z))  
  
  for(i in 1:nsim) {
    
    #Betas
    
    betas.sim <- matrix(muproposal(Y,X,Z,betas.ind,gammas.ind,bpri,Bpri),nrow=ncol(X))
    gammas.sim <- matrix(gammaproposal(Y,X,Z,betas.ind,gammas.ind,gpri,Gpri),nrow=ncol(Z))
    
    q1.mu <- mukernel(X,Z,Y,betas.sim,betas.ind,gammas.sim,bpri,Bpri)
    q2.mu <- mukernel(X,Z,Y,betas.ind,betas.sim,gammas.ind,bpri,Bpri)
    p1.mu<-dpostb(X,Z,Y,betas.sim,gammas.ind,bpri,Bpri)
    p2.mu<-dpostb(X,Z,Y,betas.ind,gammas.ind,bpri,Bpri)
    
    q1.gamma <- gammakernel(X,Z,Y,gammas.sim,betas.sim,gammas.ind,gpri,Gpri)
    q2.gamma <- gammakernel(X,Z,Y,gammas.ind,betas.ind,gammas.sim,gpri,Gpri)
    p1.gamma<-dpostg(X,Z,Y,betas.ind,gammas.sim,gpri,Gpri)
    p2.gamma<-dpostg(X,Z,Y,betas.ind,gammas.ind,gpri,Gpri)  
    
    Mu.val<-min(1,((p1.mu/p2.mu)*(q1.mu/q2.mu)))
    u<-runif(1)
    if (u <=Mu.val) {
      betas.ind <- betas.sim
      ind1[i] = 1
    }
    
    beta.mcmc[i,]<-betas.ind
    beta.mcmc <- as.ts(beta.mcmc)
    
    Gamma.val<-min(1,((p1.gamma/p2.gamma)*(q1.gamma/q2.gamma)))
    u<-runif(1)
    if (u <=Gamma.val) {
      gammas.ind <- gammas.sim
      ind2[i] = 1
    }
    gamma.mcmc[i,]<-gammas.ind
    gamma.mcmc <- as.ts(gamma.mcmc)
    
    if (i%%1000 == 0)
      cat("Burn-in iteration : ", i, "\n")
  }
  
  
  tburn <- nsim*burn
  extr <- seq(0,(nsim-tburn),jump)
  
  betas.burn <-as.matrix(beta.mcmc[(tburn+1):nrow(beta.mcmc),])
  gammas.burn <-as.matrix(gamma.mcmc[(tburn+1):nrow(gamma.mcmc),])
  
  beta.mcmc.auto <- as.matrix(betas.burn[extr,])
  beta.mcmc.auto <- as.ts(beta.mcmc.auto)
  gamma.mcmc.auto <- as.matrix(gammas.burn[extr,])
  gamma.mcmc.auto <- as.ts(gamma.mcmc.auto)
  
  
  if (graph1==TRUE) {
    
    for(i in 1:ncol(X)){
      dev.new()
      ts.plot(beta.mcmc[,i], main=paste("Complete chain for beta",i), xlab="number of iterations", ylab=paste("parameter beta",i))
    }
    
    for(i in 1:ncol(Z)){
      dev.new()
      ts.plot(gamma.mcmc[,i], main=paste("Complete chain for gamma",i), xlab="number of iterations", ylab=paste("parameter gamma",i))
      
    }
    
  } else{
  }
  
  
  if (graph2==TRUE) {
    
    for(i in 1:ncol(X)){
      dev.new()
      ts.plot(beta.mcmc.auto[,i], main=paste("Burn chain for beta",i), xlab="number of iterations", ylab=paste("parameter beta",i))
      
    }
    
    for(i in 1:ncol(Z)){
      dev.new()
      ts.plot(gamma.mcmc.auto[,i], main=paste("Burn chain for gamma",i), xlab="number of iterations", ylab=paste("parameter gamma",i))
      
    }
    
  } else{
  }
  
  #Beta y Gamma estimations
  Bestimado <- colMeans(beta.mcmc.auto)
  Gammaest <- colMeans(gamma.mcmc.auto)
  
  #estandar errors of beta and gamma
  DesvBeta <- matrix(apply(beta.mcmc.auto,2,sd))
  DesvGamma <- matrix(apply(gamma.mcmc.auto,2,sd))
  
  #Precision
  phi <- exp(Z%*%Gammaest)
  
  #estimate values of the dependent variable
  yestimado = exp(X%*%Bestimado)/(1+exp(X%*%Bestimado))
  
  #estimate variance of the dependent variable
  variance<-(yestimado*(1-yestimado))/(phi+1)
  
  #Residuals 
  residuales = as.matrix(Y) - yestimado
  
  B1 <- matrix(0, ncol(X),1)
  B2 <- matrix(0, ncol(X),1)
  
  
  #Credibility intervals for beta
  for(i in 1:ncol(X)){
    B1[i,]<-quantile(beta.mcmc.auto[,i],0.025)
    B2[i,]<-quantile(beta.mcmc.auto[,i],0.975)
    B <- cbind(B1,B2)
  }
  
  # Credibility intervals for gamma
  
  G1 <- matrix(0, ncol(Z),1)
  G2 <- matrix(0, ncol(Z),1)
  
  for(i in 1:ncol(Z)){
    G1[i,]<-quantile(gamma.mcmc.auto[,i],0.025)
    G2[i,]<-quantile(gamma.mcmc.auto[,i],0.975)
    G <- cbind(G1,G2)
  }
   list(Bestimado=Bestimado,Gammaest=Gammaest,X=X,Z=Z,DesvBeta=DesvBeta, DesvGamma=DesvGamma, B=B, G=G, yestimado=yestimado, residuales=residuales, phi=phi, variance=variance, beta.mcmc=beta.mcmc, gamma.mcmc=gamma.mcmc, beta.mcmc.auto=beta.mcmc.auto, gamma.mcmc.auto=gamma.mcmc.auto, Y = Y,X=X)
}
