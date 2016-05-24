.familyY.gaussian <- function (familyY,k,data,Y,n,z,formulaY,...){
  beta <- VarFunY <- VarY <- muY <- PY <- dispY    <- NULL
  lmodelY <- list()
  for(h in 1:k){
    modelY <- do.call(glm, 
    list(formula=formulaY,data=data,family=familyY,weights=z[,h],y=FALSE,model=FALSE,x=FALSE))
    sigma.hat <-sqrt(crossprod(sqrt(z[,h]/sum(z[,h]))*(Y-fitted(modelY))))
    beta      <- rbind(beta,coef(modelY))
    muY       <- cbind(muY,fitted(modelY)) 
    dispY     <- rbind(dispY,sigma.hat^2)
    VarFunY   <- cbind(VarFunY,rep(1,n))
    VarY      <- cbind(VarY,sigma.hat^2*rep(1,n))
    PY        <- cbind(PY,dnorm(Y,mean=fitted(modelY),sd=sigma.hat))              
    lmodelY[[h]] <- modelY
  }
  dimnames(VarFunY) <- dimnames(VarY) <- dimnames(muY) <- dimnames(PY)  <- list(1:n,paste("comp.",1:k,sep=""))
  return(list(lmodelY=lmodelY,muY=muY,dispY=as.vector(dispY),VarFunY=VarFunY,VarY=VarY,PY=PY))
}
.familyY.poisson <- function (familyY,k,data,Y,n,z,formulaY,...){
  beta <- VarFunY <- VarY <- muY <- PY <- dispY    <- NULL
  lmodelY <- list()
  for(h in 1:k){
    modelY <- do.call(glm, 
                      list(formula=formulaY,data=data,family=familyY,weights=z[,h],y=FALSE,model=FALSE,x=FALSE,
                           control=list(trace=FALSE,epsilon=1e-14)))
    beta      <- rbind(beta,coef(modelY))
    muY       <- cbind(muY,fitted(modelY)) 
    dispY     <- rbind(dispY,1)
    VarFunY   <- cbind(VarFunY,fitted(modelY))
    VarY      <- cbind(VarY,fitted(modelY))
    PY        <- cbind(PY,dpois(Y,lambda=fitted(modelY)))              
    lmodelY[[h]] <- modelY
  }
  dimnames(VarFunY) <- dimnames(VarY) <- dimnames(muY) <- dimnames(PY)  <- list(1:n,paste("comp.",1:k,sep=""))
  return(list(lmodelY=lmodelY,dispY=as.vector(dispY),VarFunY=VarFunY,VarY=VarY,PY=PY))
}
.familyY.binomial <- function (familyY,k,data,Y,n,z,mY,formulaY,...){
  # The response for a binomial distribution can be specified in a couple of ways.
  # If there is only one trial then the response can be given as a factor, with the first level
  # being "success", and all the others failure. So mY <- 1
  # The alternative is two columns (cbind'ed together) with successes in the first, 
  # and failures in the second.
 
  lmodelY <- list()
  beta <- NULL
  dispY <- numeric(k)
  VarFunY <- VarY <- muY <- PY  <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  for(h in 1:k){
    modelY <- do.call(glm, 
    list(formula=formulaY,data=data,family=familyY,weights=z[,h],y=FALSE,model=FALSE,x=FALSE,
         control=list(trace=FALSE,epsilon=1e-14)))
    
    beta      <- rbind(beta,coef(modelY))
    muY[,h]   <- mY*fitted(modelY)
    dispY[h]  <- 1
    VarFunY[,h]   <- mY*fitted(modelY)*(1-fitted(modelY))
    VarY[,h]      <- mY*fitted(modelY)*(1-fitted(modelY))
    PY[,h]        <- dbinom(Y[,1],size=mY,prob=fitted(modelY))
    lmodelY[[h]] <- modelY
  }
  return(list(lmodelY=lmodelY,dispY=as.vector(dispY),VarFunY=VarFunY,VarY=VarY,PY=PY))
}
.familyY.student.t <- function (familyY,k,data,Y,n,z,vY,t_df,formulaY,...){
  beta <- VarFunY <- VarY <- muY <- PY <- dispY <-   sig  <- NULL
  lmodelY <- list()
  zvY <- z*vY
  for(h in 1:k){
    modelY <- do.call(glm, 
                      list(formula=formulaY,data=data,family=familyY,weights=zvY[,h],y=FALSE,model=FALSE,x=FALSE,
                           control=list(trace=FALSE,epsilon=1e-14)))
    t_df[h]    <- finddf(z[,h],vY[,h],dfold=t_df[h]) #,mindf=mindf,maxdf=maxdf
    sig       <- cbind(sig,crossprod(sqrt(zvY[,h]/sum(z[,h]))*(Y-fitted(modelY))))
    
    beta      <- rbind(beta,coef(modelY))
    muY       <- cbind(muY,fitted(modelY)) 
    dispY     <- rbind(dispY,t_df[h]/(t_df[h]-2)*sig[h])
    VarFunY   <- cbind(VarFunY,rep(1,n))
    VarY      <- cbind(VarY,dispY[h]*VarFunY[,h])
    PY        <- cbind(PY,denst(Y,mu=fitted(modelY),sd=sig[h],df=t_df[h]))
    lmodelY[[h]] <- modelY
  }
  dimnames(VarFunY) <- dimnames(VarY) <- dimnames(muY) <- dimnames(PY)  <- list(1:n,paste("comp.",1:k,sep=""))
  
  return(list(lmodelY=lmodelY,dispY=as.vector(dispY),VarFunY=VarFunY,VarY=VarY,PY=PY, sig=sig, t_df=t_df,zvY=zvY))
}
.familyY.Gamma <- function (familyY,k,data,Y,n,z,formulaY,...){
  f <- function(par,Y,modelY,weights){
    l <- -sum(weights*dgam(x=Y, mu = (fitted(modelY)), nu = par, log = TRUE))
  }
  beta <- NULL
  lmodelY <- list()
  nuY  <- dispY <- numeric(k)
  VarFunY <- VarY <- muY <- PY  <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  for(h in 1:k){
    modelY <- do.call(glm, 
                      list(formula=formulaY,data=data,family=familyY,weights=z[,h],y=FALSE,model=FALSE,x=FALSE,
                           control=list(trace=FALSE,epsilon=1e-14)))
    #nuY[h]    <- MASS::gamma.shape(modelY,verbose=TRUE)$alpha
    nuY[h] <- optimize(f=f, interval=c(0,1000),Y=Y,modelY=modelY, weights=z[,h])$minimum
    beta      <- rbind(beta,coef(modelY))
    muY[,h]   <- fitted(modelY)
    dispY[h]     <- 1/nuY[h]
    VarFunY[,h]   <- muY[,h]^2
    VarY[,h]      <- dispY[h]*VarFunY[,h]
    PY[,h]        <- dgam(Y, mu = muY[,h], nu = nuY[h])
    lmodelY[[h]] <- modelY
  }
  return(list(lmodelY=lmodelY,dispY=as.vector(dispY),VarFunY=VarFunY,VarY=VarY,nuY=nuY,PY=PY))
} 
.familyY.inverse.gaussian <- function (familyY,k,data,Y,n,z,formulaY,...){
  f <- function(par,Y,modelY,weights){
    l <- -sum(weights*dig(x=Y, mu=fitted(modelY), var=par, log=TRUE))
  }
  beta <- NULL
  lmodelY <- list()  
  dispY <- numeric(k)
  VarFunY <- VarY <- muY <- PY  <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  for(h in 1:k){
    modelY <- do.call(glm, 
                      list(formula=formulaY,data=data,family=familyY,weights=z[,h],y=FALSE,model=FALSE,x=FALSE,
                           control=list(trace=FALSE,epsilon=1e-14)))
    dispY[h] <- optimize(f=f, interval=c(0,1000), Y=Y, modelY=modelY, weights=z[,h])$minimum
    beta     <- rbind(beta,coef(modelY))
    muY[,h]  <- fitted(modelY)
    VarFunY[,h]   <- muY[,h]^3
    VarY[,h]      <- dispY[h]*VarFunY[,h]
    PY[,h]        <- dig(Y, mu = muY[,h], var = dispY[h])
    lmodelY[[h]] <- modelY
  }
  return(list(lmodelY=lmodelY,dispY=as.vector(dispY),VarFunY=VarFunY,VarY=VarY,PY=PY))
} 