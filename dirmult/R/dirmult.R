## Help-functions for computing the score-function and log-likelihood function
uu <- function(x,t) 1/(t+x-1)
ff <- function(x,t) 1/(t+x-1)^2
ll <- function(x,t) log(t+x-1)

## Computes Beta-Binomial probabilities
dbbin.ab <- function(x,n,a,b){
  res <- lchoose(n,x)
  if(x==0) res <- res else res <- res+sum(unlist(lapply(list(1:x),ll,t=a)))
  if(x==n) res <- res else res <- res+sum(unlist(lapply(list(1:(n-x)),ll,t=b)))
  res <- res-sum(unlist(lapply(list(1:n),ll,t=a+b)))
  exp(res)
}

## Score-function
u <- function(x,t){
  nc <- ncol(x)
  S <- rep(0,nc)
  ts <- sum(t)
  for(j in 1:nrow(x)){
    Sn <- sum(unlist(lapply(list(1:(rowSums(x)[j])),uu,t=ts)))
    for(i in 1:nc){
      if(x[j,i]==0) Sij <- 0
      else Sij <- sum(unlist(lapply(list(1:(x[j,i])),uu,t=t[i])))
      S[i] <- S[i] - Sn + Sij
    }
  }
  S
}

## Score-function when Lagrange multiplier is envoked
profU <- function(x,t,tp){
  K <- length(t)
  S <- rep(0,K+1)
  nc <- ncol(x)
  ts <- sum(t)
  for(j in 1:nrow(x)){
    Sn <- sum(unlist(lapply(list(1:(rowSums(x)[j])),uu,t=ts)))
    for(i in 1:nc){
      if(x[j,i]==0) Sij <- 0
      else Sij <- sum(unlist(lapply(list(1:(x[j,i])),uu,t=t[i])))
      S[i] <- S[i] - Sn + Sij
    }
  }
  S[1:K] <- S[1:K]-tp
  S[K+1] <- tp-sum(t)
  S
}

## Score-function with several Lagrange multiplier for testing equal thetas
equalU <- function(x,t,tp,l){
  K <- length(t)
  S <- rep(0,K+1)
  nc <- ncol(x)
  ts <- sum(t)
  for(j in 1:nrow(x)){
    Sn <- sum(unlist(lapply(list(1:(rowSums(x)[j])),uu,t=ts)))
    for(i in 1:nc){
      if(x[j,i]==0) Sij <- 0
      else Sij <- sum(unlist(lapply(list(1:(x[j,i])),uu,t=t[i])))
      S[i] <- S[i] - Sn + Sij
    }
  }
  S[1:K] <- S[1:K]-l
  S[K+1] <- tp-sum(t)
  S
}

## Computing the observed Fisher Information Matrix
obsfim <- function(x,t){
  nc <- ncol(x)
  D <- rep(0,nc)
  od <- 0
  ts <- sum(t)
  for(j in 1:nrow(x)){
    od <- od + sum(unlist(lapply(list(1:(rowSums(x)[j])),ff,t=ts)))
    for(i in 1:nc){
      if(x[j,i]==0) Dij <- 0
      else Dij <- sum(unlist(lapply(list(1:(x[j,i])),ff,t=t[i])))
      D[i] <- D[i] + Dij
    }
  }
  F <- matrix(od,nc,nc)
  diag(F) <- diag(F)-D
  F
}

## Computing the expected Fisher Information Matrix
expfim <- function(x,t){
  Sn <- rowSums(x)
  J <- nrow(x)
  K <- ncol(x)
  inner <- matrix(0,1,K)
  od <- 0
  ts <- sum(t)
  for(j in 1:J){
    P <- matrix(0,Sn[j],K)
    R <- matrix(0,Sn[j],K)
    for(r in Sn[j]:1){
      for(k in 1:K){
        if(r==Sn[j]) P[r,k] <- dbbin.ab(r,Sn[j],t[k],ts-t[k])
        else P[r,k] <- P[r+1,k] + dbbin.ab(r,Sn[j],t[k],ts-t[k])
        R[r,k] <- 1/(t[k]+r-1)^2
      }
    }
    inner <- inner + colSums(P*R)
    od <- od + sum(unlist(lapply(list(1:(Sn[j])),ff,t=ts)))
  }
  F <- matrix(-1*od,K,K)
  diag(F) <- diag(F)+inner
  F
}

## Computes the log-likelihood
loglik <- function(x,t){
  l <- 0
  ts <- sum(t)
  nc <- ncol(x)
  for(j in 1:nrow(x)){
    l <- l - sum(unlist(lapply(list(1:(rowSums(x)[j])),ll,t=ts)))
# maybe not    lij <- 0 ## NEW LINE: This initializes lij to zero for each row. Otherwise lij is a cummulant. ##
    for(i in 1:nc){
      if(x[j,i]==0) lij <- 0
      else lij <- sum(unlist(lapply(list(1:(x[j,i])),ll,t=t[i])))
      l <- l + lij
    }
  }
  l
}

## Computes: I(pi,theta) = t(D)%*%I(gamma)%*%D, where D is the Jacobi matrix
thetafim <- function(t,f){
  K <- length(t)
  D <- matrix(0,K,K)
  pi <- t/sum(t)
  theta <- 1/(sum(t)+1)
  diag(D) <- (1-theta)/theta
  D[K,] <- -(1-theta)/theta
  D[,K] <- -1*pi/(theta^2) 
  D[K,K] <- D[K,K] 
  t(D)%*%f%*%D
}

## Estimate parameters in the Dirichlet-Multinomial distribution
dirmult <- function(data,init,initscalar,epsilon=10^(-4),trace=TRUE,mode){
  data <- data[rowSums(data)!=0,colSums(data)!=0]
  if(missing(initscalar)){
    mom <- weirMoM(data)
    if(mom<=0) mom <- 0.005
    initscalar <- (1-mom)/mom
  }
  if(missing(init)) gamma <- colSums(data)/sum(data)*initscalar
  else gamma <- init
  if(missing(mode)) mode <- "obs"
  if(!is.element(mode,c("obs","exp"))){
    message(paste("Warning: Mode '",mode,"' not valid\n",sep=""))
    mode <- "obs"
  }
  lik1 <- 0
  lik2 <- epsilon*10
  ite <- 1
  gamite <- 0
  conv <- TRUE
  # Iterations
  while(conv){
    if(abs(lik2-lik1)<epsilon) conv <- FALSE
    if(mode=="exp") fim <- expfim(data,gamma)
    else if(mode=="obs") fim <- -1*obsfim(data,gamma)
    lik1 <- loglik(data,gamma)
    # Updates parameter estimates
    gamma <- gamma + solve(fim)%*%u(data,gamma)
    gamma[gamma<0] <- 0.01 # Negative gamma_j are set to 0.01
    if(any(gamma<0) & (gamite%%10)==0){ print(gamma); gamite <- gamite+1}
    if(trace) message(paste("Iteration ",ite,": Log-likelihood value: ",lik1,sep=""))
    gams <- paste(" Gamma",1:length(gamma),sep="")
    lik2 <- loglik(data,gamma)
    ite <- ite+1
  }
  sumgam <- sum(gamma)
  theta <- 1/(sumgam+1)
  pi <- as.numeric(gamma/sumgam)
  names(pi) <- dimnames(data)[[2]]
  list(loglik=lik1,ite=ite-1,gamma=as.numeric(gamma),pi=pi,theta=theta)
}

## Generating a summary table with estimates and std. errors for MLE and MoM
dirmult.summary <- function(data,fit,expectedFIM=FALSE){
  K <- ncol(data)
  J <- nrow(data)
  if(expectedFIM) fim <- expfim(data,fit$gamma)
  else fim <- -obsfim(data,fit$gamma)
  FIM <- thetafim(fit$gamma,fim)
  invFIM <- solve(FIM)
  stdMLE <- sqrt(diag(invFIM))
  FIMpi <- invFIM[-K,-K]
  stdPik <- sqrt(sum(stdMLE[-K]^2)+2*sum(FIMpi[upper.tri(FIMpi)]))
  MoM <- colSums(data)/sum(data)
  Sn <- rowSums(data)
  MSP <- (J-1)^(-1)*sum(rowSums((data/rowSums(data)-matrix(rep(MoM,J),J,K,byrow=T))^2)*Sn)
  MSG <- (sum(data)-J)^(-1)*sum(rowSums(data/rowSums(data)*(1-data/rowSums(data)))*Sn)
  nc <- 1/(J-1)*(sum(Sn)-sum(Sn^2)/sum(Sn))
  weir <- weirMoM(data,se=TRUE)
  MoM.wh <- weir$theta
  MoM.se <- sqrt(colSums((data/rowSums(data)-matrix(rep(MoM,J),J,K,byrow=T))^2)/(J-1))
  MoM.se <- c(MoM.se,weir$se)
  res <- data.frame("Variable"=c(paste("pi",dimnames(data)[[2]],sep=":"),"Theta"),
                    "MLE"=c(fit$pi,fit$theta),
                    "se.MLE"=c(stdMLE[-K],stdPik,stdMLE[K]),
                    "MoM"=c(MoM,MoM.wh), 
                    "se.MOM"=MoM.se)
  res[order(res$Variable),]
}

## Computes the MoM estimate of theta (and std. error)
weirMoM <- function(data,se=FALSE){
  K <- ncol(data)
  J <- nrow(data)
  MoM <- colSums(data)/sum(data)
  Sn <- rowSums(data)
  MSP <- (J-1)^(-1)*sum(rowSums((data/rowSums(data)-matrix(rep(MoM,J),J,K,byrow=T))^2)*Sn)
  MSG <- (sum(data)-J)^(-1)*sum(rowSums(data/rowSums(data)*(1-data/rowSums(data)))*Sn)
  nc <- 1/(J-1)*(sum(Sn)-sum(Sn^2)/sum(Sn))
  MoM.wh <- (MSP-MSG)/(MSP+(nc-1)*MSG)
  if(se){
    ## Formula by Li, ref in Weir-Hill 2002
    std.er <- sqrt(2*(1-MoM.wh)^2/(J-1)*((1+(nc-1)*MoM.wh)/nc)^2)
    list(theta=MoM.wh,se=std.er)
    }
  else MoM.wh
}

## Estimates the value of the profile log-likelihood in value 'theta'
estProfLogLik <- function(data,theta,epsilon=10^(-4),trace=TRUE,initPi,maxit=1000){
  gamplus <- (1-theta)/theta
  data <- data[,colSums(data)!=0]
  K <- ncol(data)
  if(!missing(initPi)) gamma <- initPi*gamplus
  else gamma <- colSums(data)/sum(data)*gamplus
  gamlambda <- c(gamma,1)
  lik1 <- 0
  lik2 <- epsilon*10
  ## iterations:
  ite <- 1
  gamite <- 0
  conv <- TRUE
  while(conv){
    if(abs(lik2-lik1)<epsilon) conv <- FALSE
    if(ite>maxit) return(NULL)
    fimGam <- obsfim(data,gamma)
    fim <- matrix(-1,K+1,K+1)
    fim[K+1,K+1] <- 0
    fim[1:K,1:K] <- fimGam
    lik1 <- loglik(data,gamma)+gamlambda[K+1]*(gamplus-sum(gamma))
    # Updates parameter estimates
    gamlambda <- gamlambda - solve(fim)%*%profU(data,gamma,gamplus)
    gamma <- gamlambda[1:K]
    if(any(gamma<0)){
      if((gamite%%10)==0) print(gamma)
      gamite <- gamite+1
    }
    gamma[gamma<0] <- 0.001 # Negative gamma_j are set to 0.001
    if(trace) message(paste("Iteration ",ite,": Log-likelihood value: ",lik1,sep=""))
    if(ite>50) message(paste("Iteration ",ite,": Log-likelihood value: ",lik1,sep=""))
    gams <- paste(" Gamma",1:length(gamma),sep="")
    lik2 <- loglik(data,gamma)+gamlambda[K+1]*(gamplus-sum(gamma))
    ite <- ite+1
  }
  gamma <- gamlambda[1:K]
  sumgam <- sum(gamma)
  pi <- as.numeric(gamma/sumgam)
  names(pi) <- dimnames(data)[[2]]
  list(loglik=lik1,ite=ite-1,gamma=as.numeric(gamma),pi=pi,theta=1/(sumgam+1),lambda=gamlambda[K+1])
}

## Estimates the profile log-likelihood for a defined grid of values
gridProf <- function(data,theta,from,to,len){
  step <- theta+seq(from=from,to=to,len=len)
  res <- data.frame(theta=step,loglik=rep(0,len))
  for(i in 1:len) res$loglik[i] <- estProfLogLik(data,step[i],trace=FALSE)$loglik
  res
}

## Estimates the profile log-likelihood such that difference in loglik is max delta
adapGridProf <- function(data,delta,stepsize=50){
  mle <- dirmult(data,trace=FALSE,epsilon=10^(-8))
  if(is.na(mle$theta)) print("MLE theta=NA")
  step <- mle$theta/stepsize ## use /100 or /1000 for better precision in CI
  res <- data.frame(theta=mle$theta,loglik=mle$loglik)
  stopp <- TRUE
  k <- 1
  pip <- mle$pi
  pim <- mle$pi
  # Algorithm works its way to the right of MLE
  while(stopp){
    tmpp <- mle$theta+k*step
    llp <- estProfLogLik(data,tmpp,trace=FALSE,initPi=pip)
    if(is.null(llp)) return(NULL)
    if(llp$ite>300) return(NULL)
    pip <- llp$pi
    llp <- llp$loglik
    res <- rbind(res,c(tmpp,llp))
    if(is.na(llp) | abs(llp-mle$loglik)>delta) stopp <- FALSE
    else k <- k+1
  }
  stopm <- TRUE
  l <- 1
  # Algorithm works its way to the left of MLE
  while(stopm){
    tmpm <- mle$theta-l*step
    llm <- estProfLogLik(data,tmpm,trace=FALSE,initPi=pim)
    if(is.null(llm)) return(NULL)
    if(llm$ite>300) return(NULL)
    pim <- llm$pi
    llm <- llm$loglik
    res <- rbind(res,c(tmpm,llm))
    if(is.na(llm) | abs(llm-mle$loglik)>delta) stopm <- FALSE
    else l <- l+1
  }
  res[order(res$theta),]
}

## Computes the log-likelihood function assuming equal theta for all tables in the list 'data' 
equalTheta <- function(data,theta,epsilon=10^(-4),trace=TRUE,initPi,maxit=1000){
  gamplus <- (1-theta)/theta
  data <- lapply(data,function(x) x[,colSums(x)!=0])
  L <- length(data)
  K <- unlist(lapply(data,ncol))
  KK <- K+1 ## dim(gamma) + dim(lambda)
  sKK <- sum(KK)
  cKK <- c(0,cumsum(KK))
  if(!missing(initPi)) gamma <- lapply(initPi,function(x,t) x*t,t=gamplus)
  else gamma <- lapply(data,function(x) colSums(x)/sum(x)*gamplus)
  gamlambda <- lapply(gamma,function(x) c(x,1))
  fimGam <- as.list(rep(0,L))
  invfimGam <- as.list(rep(0,L))
  fim <- as.list(rep(0,L))
  invfim <- as.list(rep(0,L))
  deninvfim <- as.list(rep(0,L))
  numinvfim <- as.list(rep(0,L))
  scorevector <- as.list(rep(0,L))
  lik1 <- 0
  lik2 <- epsilon*10
  ite <- 1
  conv <- TRUE
  while(conv){
    if(abs(lik2-lik1)<epsilon) conv <- FALSE
    if(ite>maxit) return(NULL)    
    lik1 <- 0
    FIM <- matrix(0,sKK+1,sKK+1)
    FIMtest <- matrix(0,sKK+1,sKK+1)
    lambda <- unlist(lapply(gamlambda,function(x) x[length(x)]))
    for(l in 1:L){
      # Fits for each data table in list 'data'
      fimGam[[l]] <- obsfim(data[[l]],gamma[[l]])
      invfimGam[[l]] <- solve(fimGam[[l]])
      fim[[l]] <- matrix(-1,KK[l],KK[l])
      fim[[l]][KK[l],KK[l]] <- 0
      fim[[l]][1:K[l],1:K[l]] <- fimGam[[l]]
      invfim[[l]] <- solve(fim[[l]])
      deninvfim[[l]] <- (-1)/(matrix(1,1,K[l])%*%invfimGam[[l]]%*%matrix(1,K[l],1))
      numinvfim[[l]] <- invfim[[l]]%*%matrix(c(rep(0,K[l]),1),K[l]+1,1)
      FIM[(cKK[l]+1):cKK[l+1],(cKK[l]+1):cKK[l+1]] <- invfim[[l]]
      FIMtest[(cKK[l]+1):cKK[l+1],(cKK[l]+1):cKK[l+1]] <- fim[[l]]
      scorevector[[l]] <- equalU(data[[l]],gamma[[l]],gamplus,lambda[l])
      lik1 <- lik1 + loglik(data[[l]],gamma[[l]])+lambda[l]*(gamplus-sum(gamma[[l]]))
    }
    # Computes overall parameters and log-likelihood
    kk <- rep(1,length(K)*2)
    kk[rep(c(T,F),L)] <- K
    FIMtest[sKK+1,1:sKK] <- rep(rep(c(0,1),L),kk)
    FIMtest[1:sKK,sKK+1] <- rep(rep(c(0,1),L),kk)
    FIMtest[sKK+1,sKK+1] <- 0
    Ax <- unlist(numinvfim)
    xAx <- sum(unlist(deninvfim))
    FIM[1:sKK,1:sKK] <- FIM[1:sKK,1:sKK]-Ax%*%t(Ax)/xAx
    FIM[1:sKK,sKK+1] <- Ax/xAx
    FIM[sKK+1,1:sKK] <- t(Ax)/xAx
    FIM[sKK+1,sKK+1] <- (-1)/xAx
    gammalambda <- c(unlist(gamlambda),gamplus)
    uvector <- c(unlist(scorevector),sum(lambda))
    gammalambda <- gammalambda - FIM%*%uvector
    gamplus <- gammalambda[sKK+1]
    gamlambda <- split(gammalambda[-(sKK+1)],as.factor(rep(1:L,KK)))
    gamma <- lapply(gamlambda,function(x) x[-length(x)])
    if(any(unlist(gamma)<0)){ ## NEGATIVE ENTRY IN GAMMA VECTOR ##
      neg <- (1:L)[unlist(lapply(gamma,function(x) any(x<0)))]
      print(gamma[neg])
      gamma <- lapply(gamma,function(x){ x[x<0] <- 0.01; x}) # set negative entries to 0.01
    }
    if(trace) message(paste("Iteration ",ite,": Log-likelihood value: ",lik1,sep=""))
    lik2 <- 0
    for(l in 1:L) lik2 <- lik2 + loglik(data[[l]],gamma[[l]])+lambda[l]*(gamplus-sum(gamma[[l]]))
    ite <- ite+1
  }
  pi <- lapply(gamma,function(x) x/sum(x))
  for(l in 1:L) names(pi[[l]]) <- dimnames(data[[l]])[[2]]
  theta <- lapply(gamma,function(x) 1/(1+sum(x)))
  list(loglik=lik1,ite=ite-1,gamma=gamma,pi=pi,theta=theta,lambda=lambda)
}

#log-like of multinomial
mnloglik <- function(x){
  x <- x[,colSums(x)!=0]
  p <- colSums(x)/sum(x)
  sum(x*rep(log(p),each=nrow(x)))
}

## Simulates under H_0: theta=0 
nullTest <- function(data,m=1000,prec=6){
  pi.null <- colSums(data)/sum(data)
  dats <- replicate(m,data,simplify=FALSE)
  res <- data.frame(mle=rep(0,m+1),dm=rep(0,m+1),mom=rep(0,m+1),mn=rep(0,m+1))
  rs <- rowSums(data)
  nr <- nrow(data)
  for(i in 1:m){
    for(j in 1:nr) dats[[i]][j,] <- rmultinom(1,rs[j],pi.null)
    tmp <- unlist(dirmult(dats[[i]],trace=FALSE,epsilon=10^(-prec))[c("theta","loglik")])
    res[i,] <- c(tmp,weirMoM(dats[[i]]),mnloglik(dats[[i]]))
  }
  tmp <- unlist(dirmult(data,trace=FALSE,epsilon=10^(-prec))[c("theta","loglik")])
  res[m+1,] <- c(tmp,weirMoM(data),mnloglik(data))
  list(data=dats,res=res)
}

rdirichlet <- function(n=1,alpha){
  Gam <- matrix(0,n,length(alpha))
  for(i in 1:length(alpha)) Gam[,i] <- rgamma(n,shape=alpha[i])
  Gam/rowSums(Gam)
}

simPop <- function(J=10,K=20,n,pi,theta){
  if(length(n)==1) n <- rep(n,J)
  if(missing(pi)) pi <- rnorm(K,mean=14,sd=4)
  else K <- length(pi)
  pi <- pi/sum(pi)
  P <- rdirichlet(J,pi*(1-theta)/theta)
  X <- matrix(0,J,K)
  for(i in 1:J) X[i,] <- rmultinom(1,n[i],P[i,])
  list(theta=theta,pi=pi,data=X)
}

