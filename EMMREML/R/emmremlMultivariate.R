
emmremlMultivariate<-function(Y,X, Z, K,varBhat=FALSE,varGhat=FALSE, PEVGhat=FALSE, test=FALSE,tolpar=1e-06, tolparinv=1e-06){
	
   Z<-t(Z)
 ECM1<-function(ytl, xtl, Vgt,Vet,Bt, deltal){
 Vlt=deltal*Vgt+Vet
 	invVlt<-solve(Vlt+tolparinv*diag(d))
 return(list(Vlt=Vlt, gtl=deltal*Vgt%*%invVlt%*%(ytl-Bt%*%xtl), Sigmalt=deltal*Vgt-deltal*Vgt%*%invVlt%*%(deltal*Vgt)))

}

wrapperECM1<-function(l){
ytl<-Yt[,l]
xtl<-Xt[,l]
 deltal<-eigZKZt$values[l]

return( ECM1(ytl=ytl, xtl=xtl, Vgt=Vgt,Vet=Vet,Bt=Bt, deltal=deltal))
 }



Vgfunc<-function(l){
	Vgl<-tcrossprod(outfromECM1[[l]]$gtl)
		return((1/n)*(1/eigZKZt$values[l])*(Vgl + outfromECM1[[l]]$Sigmalt))
		}
		
Vefunc<-function(l){
		etl <- Yt[,l] - Bt%*%Xt[,l] - outfromECM1[[l]]$gtl
		return((1/n)*((tcrossprod(etl)+ outfromECM1[[l]]$Sigmalt)))
		}



 if (sum(is.na(Y))==0){
  N<-nrow(K)
  KZt<-tcrossprod(K,Z)
  ZKZt<-Z%*%KZt
   
 eigZKZt = eigen(ZKZt)
 n<-nrow(ZKZt)
  d<-nrow(Y)
 Yt = Y%*%eigZKZt$vectors
 Xt = X%*%eigZKZt$vectors

Vgt =cov(t(Y))/2
 Vet =cov(t(Y))/2
  XttinvXtXtt<-t(Xt)%*%solve(tcrossprod(Xt))
Bt<-Yt%*%XttinvXtXtt
Vetm1<-Vet
repeat{
outfromECM1<-lapply(1:n, wrapperECM1)
Vetm1<-Vet

	Gt=sapply(outfromECM1, function(x) {cbind(x$gtl)})
	Bt = (Yt - Gt) %*% XttinvXtXtt
	
	listVgts <- lapply(1:n,Vgfunc)
	Vgt<-Reduce('+', listVgts)
	listVets <- lapply(1:n,Vefunc)
	Vet<-Reduce('+', listVets)
	convnum<-abs(sum(diag(Vet - Vetm1)))/abs(sum(diag(Vetm1)))
	convcond<-tryCatch({convnum<tolpar}, error=function(e){return(FALSE)})
	if(convcond){break}
}


HobsInv<-solve(kronecker(ZKZt,Vgt)+kronecker(diag(n),Vet)+tolparinv*diag(d*n))
ehat<-matrix(Y - Bt%*%X,ncol=1, byrow=F)
HobsInve<-HobsInv%*%ehat
varvecG<-kronecker(K,Vgt)
gpred<-varvecG%*%(kronecker(t(Z),diag(d)))%*%HobsInve
Gpred<-matrix(gpred, nrow=nrow(Y), byrow=F)
colnames(Gpred)<-rownames(K)

  		Xforvec<-(kronecker(t(X),diag(d)))
  		Zforvec<-(kronecker((Z),diag(d)))
  		ZKforvec<-Zforvec%*%varvecG
  ####VAR U
  if (varGhat){
  	
  			P<-HobsInv-HobsInv%*%Xforvec%*%solve(crossprod(Xforvec,HobsInv%*%Xforvec), crossprod(Xforvec,HobsInv))
  			varGhat<-crossprod(ZKforvec,P)%*%ZKforvec
}
  
   if (PEVGhat){
  			
  			if (!exists("P")){P<-HobsInv-HobsInv%*%Xforvec%*%solve(crossprod(Xforvec,HobsInv%*%Xforvec), crossprod(Xforvec,HobsInv))}
  			PEVGhat<-varvecG-varGhat
}

   #varbeta
  
   if (varBhat){
  	varBhat<-solve(crossprod(Xforvec,HobsInv%*%Xforvec))
  }
  if (test){
  	XsqtestG<-matrix(Gpred,ncol=1, byrow=F)^2/diag(varGhat)
  	XsqtestG<-matrix(XsqtestG, ncol=nrow(Y), byrow=T)
  	XsqtestG<-matrix(XsqtestG,ncol=1, byrow=T)
	pGhat<-pchisq(XsqtestG,df=1, lower.tail=F,log.p=F)
	p.adjust.M <- p.adjust.methods
	p.adjGhat    <- sapply(p.adjust.M, function(meth) p.adjust(pGhat, meth))
	XsqtestB<-matrix(Bt,ncol=1, byrow=F)^2/diag(varBhat)
	pBhat<-pchisq(XsqtestB,df=1, lower.tail=F,log.p=F)
	p.adjBhat    <- sapply(p.adjust.M, function(meth) p.adjust(pBhat, meth))
  }
  if (!exists("XsqtestB")){XsqtestB<-c()}
  if (!exists("p.adjBhat")){p.adjBhat<-c()}
  if (!exists("XsqtestG")){XsqtestG<-c()}
  if (!exists("p.adjGhat")){p.adjGhat<-c()}
  if (!exists("varGhat")){varGhat<-c()}
  if (!exists("varBhat")){varBhat<-c()}
  if (!exists("PEVGhat")){PEVGhat<-c()}


return(list(Bhat=Bt,Vg=Vgt,Ve=Vet, Gpred=Gpred, XsqtestB=XsqtestB,pvalB=p.adjBhat,XsqtestG=XsqtestG,pvalG=p.adjGhat,varGhat=diag(varGhat), varBhat=diag(varBhat), PEVGhat=diag(PEVGhat)))
}
}
