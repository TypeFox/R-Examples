emmremlMultiKernel <-
function(y, X, Zlist, Klist,varbetahat=FALSE,varuhat=FALSE, PEVuhat=FALSE, test=FALSE){
	q=dim(X)[2]
  	n=length(y)
   	lz<-length(Zlist)
  
  spI<- diag(n)
  S<-spI-tcrossprod(X%*%solve(crossprod(X)),X)
  
 
 Z<-c()
 for (i in 1:lz){
  Z<-cbind(Z, Zlist[[i]])
  }
  

 
  minimfunctionouter<-function(weights=rep(1/lz, lz)){
    weights=weights/sum(weights)
    ZKZt<-matrix(0,nrow=n,ncol=n)
  for (i in 1:lz){ZKZt<-ZKZt+weights[i]*tcrossprod(Zlist[[i]]%*%Klist[[i]],Zlist[[i]])}

  offset<-log(n)
  ZKZtandoffset<-ZKZt+offset*spI
  SZKZtSandoffset<-{S%*%ZKZtandoffset}%*%S
  
  svdSZKZtSandspI<-eigen(SZKZtSandoffset, symmetric=TRUE)
  Ur<-svdSZKZtSandspI$vectors[,1:(n-q)]
  
  lambda<-svdSZKZtSandspI$values[1:(n-q)]-offset
  eta<-crossprod(Ur,y)
  minimfunc<-function(delta){(n-q)*log(sum(eta^2/{lambda+delta}))+sum(log(lambda+delta))}
  optimout<-optimize(minimfunc, lower=0, upper=10000, tol=1e-9)
  return(optimout$objective)
  }
  weights<-optim(par=rep(1/lz, lz), fn=minimfunctionouter, 
      method = "L-BFGS-B", lower = rep(0, lz), upper = rep(1,lz))$par
  Z <- c()
    for (i in 1:lz) {
        Z <- cbind(Z, Zlist[[i]])
    }

  weights<-weights/sum(weights)
  ZKZt<-matrix(0,nrow=n,ncol=n)
  Klistweighted<-Klist
  for (i in 1:lz){
  	Klistweighted[[i]] <- weights[i] * Klist[[i]]
  	ZKZt<-ZKZt+weights[i]*tcrossprod(Zlist[[i]]%*%Klist[[i]],Zlist[[i]])
  	
  	}
K <- .bdiag(Klistweighted)
    #K <- as.matrix(K)
    ZK <- as.matrix(Z %*% K)

  offset<-log(n)
  ZKZtandoffset<-ZKZt+offset*spI

  SZKZtSandoffset<-{S%*%ZKZtandoffset}%*%S
  
  svdSZKZtSandspI<-eigen(SZKZtSandoffset, symmetric=TRUE)
  Ur<-svdSZKZtSandspI$vectors[,1:(n-q)]
  
  lambda<-svdSZKZtSandspI$values[1:(n-q)]-offset
  eta<-crossprod(Ur,y)
  minimfunc<-function(delta){(n-q)*log(sum(eta^2/{lambda+delta}))+sum(log(lambda+delta))}
  optimout<-optimize(minimfunc, lower=0, upper=10000, tol=1e-9)
  deltahat<-optimout$minimum

  #optimout$minimum
  #print(deltahat)
  Hinvhat<-solve(ZKZt+deltahat*spI)
  XtHinvhat<-crossprod(X,Hinvhat)
  betahat<-solve(XtHinvhat%*%X,XtHinvhat%*%y)
  ehat<-(y-{X%*%betahat})
  Hinvhatehat<-Hinvhat%*%ehat
  
  sigmausqhat<-sum(eta^2/{lambda+deltahat})/(n-q)
  sigmaesqhat<-deltahat*sigmausqhat
  uhat<-crossprod(ZK,Hinvhatehat)
  Vinv<-(1/sigmausqhat)*Hinvhat
  namesuhat<-c()
  for (i in 1:length(Klist)){namesuhat<-c(namesuhat,paste(paste("K",i,sep="_"),colnames(Klist[[i]]),sep=""))}
  
  df <- n - q
  loglik <-  -0.5 * (optimout$objective + df + df * log(2 * pi/df))
  ####VAR U
  if (varuhat){
  			
  			P<-Vinv-Vinv%*%X%*%solve(crossprod(X,Vinv%*%X), crossprod(X,Vinv))
  			varuhat<-sigmausqhat^2*crossprod(ZK,P)%*%ZK
}
  
   if (PEVuhat){
  			
  			if (!exists("P")){P<-Vinv-Vinv%*%X%*%solve(crossprod(X,Vinv%*%X), crossprod(X,Vinv))}
  			PEVuhat<-sigmausqhat*K-varuhat
}

   #varbeta
  
   if (varbetahat){
  	varbetahat<-solve(crossprod(X,Vinv%*%X))
  }
  if (test){
  	Xsqtestu<-uhat^2/diag(varuhat)
	puhat<-pchisq(Xsqtestu,df=1, lower.tail=F,log.p=F)
	p.adjust.M <- p.adjust.methods
	p.adjuhat    <- sapply(p.adjust.M, function(meth) p.adjust(puhat, meth))
	Xsqtestbeta<-betahat^2/diag(varbetahat)
	pbetahat<-pchisq(Xsqtestbeta,df=1, lower.tail=F,log.p=F)
	p.adjbetahat    <- sapply(p.adjust.M, function(meth) p.adjust(pbetahat, meth))
  }
  if (!exists("Xsqtestbeta")){Xsqtestbeta<-c()}
  if (!exists("p.adjbetahat")){p.adjbetahat<-c()}
  if (!exists("Xsqtestu")){Xsqtestu<-c()}
  if (!exists("p.adjuhat")){p.adjuhat<-c()}
  if (!exists("varuhat")){varuhat<-c()}
  if (!exists("varbeta")){varubeta<-c()}
  if (!exists("PEVuhat")){PEVuhat<-c()}
uhat<-as.numeric(uhat)
names(uhat)<-namesuhat
  return(list(Vu=sigmausqhat,Ve=sigmaesqhat,betahat=betahat,uhat=uhat, weights=weights, Xsqtestbeta=Xsqtestbeta,pvalbeta=p.adjbetahat,Xsqtestu=Xsqtestu,pvalu=p.adjuhat,varuhat=diag(varuhat), varbetahat=diag(varbetahat), PEVuhat=diag(PEVuhat), loglik=loglik))

}
