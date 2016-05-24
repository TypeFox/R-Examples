reml=function(y,X=NULL,Z=NULL,K=NULL){
  anyNA = function(x) any(is.na(x))
  if(!is.matrix(y)){
    
    N = length(y)
    
    # Dealing with fixed effect matrix
    if(is.null(X)){X = matrix(1,N,1)}else{
      if(is.matrix(X)){if(nrow(X)!=N) stop("Fixed effect does not match dimensions of response variable")
      }else{if(class(X)=="formula"){X=model.matrix(X)}}}
    
    # Dealing with random effect
    if(is.null(K)&is.null(Z))stop("Either Z or K must be specified")
    if(is.null(K)){
      if(class(Z)=="formula"){Z=model.matrix(Z)-1}
      V=tcrossprod(Z)}
    if(is.null(Z)){V=K;Z=diag(ncol(K))}
    if(is.null(Z)!=T&&is.null(K)!=T){
      if(class(Z)=="formula"){Z=model.matrix(Z)-1}
      V=crossprod(t(Z),K);V=tcrossprod(V,Z)}
    K=V
    
    # Function starts here
    m = which(is.na(y)) # missing values
    if(any(is.na(y))){
      y=y[-m];x=X[-m,]
      k=K[m,-m];K=K[-m,-m]
      z=Z[-m,]
    }else{
      x=X;k=K;z=Z
      }
    x=as.matrix(x)
    # Defining log-REML
    loglike=function(theta){
      lambda=exp(theta)
      logdt=sum(log(lambda*delta+1))
      h=1/(lambda*delta+1)
      yy=sum(yu*h*yu);yx=matrix(0,q,1)
      xx=matrix(0,q,q);for(i in 1:q){
        yx[i]=sum(yu*h*xu[,i])
        for(j in 1:q){xx[i,j]=sum(xu[,i]*h*xu[,j])}}
      loglike = -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
      return(-loglike)}
    fixed=function(lambda){
      h=1/(lambda*delta+1)
      yy=sum(yu*h*yu)
      yx=timesVec(yu,h,xu,q)
      xx=timesMatrix(xu,h,xu,q,q)
      beta=qr.solve(xx,yx)
      sigma2=(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
      sigma2 = as.numeric(sigma2)
      var=diag((chol2inv(xx))*sigma2)
      stderr=sqrt(var)
      return(c(beta,stderr,sigma2))}
    # Eigendecomposition of K
    qq=eigen(as.matrix(K),symmetric=T)
    delta=qq[[1]]; uu=qq[[2]]; q=max(ncol(x),1)
    n=ncol(K); vp=var(y)
    yu=t(uu)%*%y; xu=t(uu)%*%x
    theta=0
    # Finding lambda through optimization
    parm=optim(par=theta,fn=loglike,method="L-BFGS-B",lower=-10,upper=10)
    lambda=exp(parm$par)
     # Variance components and fixed effect coefficient
    parmfix=fixed(lambda)
    beta=parmfix[1:q]
    sd=parmfix[(q+1):(2*q)]
    B = cbind(beta,sd)
    Ve=parmfix[2*q+1]; Vg=lambda*Ve
    h2=Vg/(Vg+Ve); VC = data.frame(Vg,Ve,h2)
    # Random effect coefficient and prediction
    # VanRaden(2008)
    if(is.matrix(x)){re =y-x%*%beta}else{re =(y-tcrossprod(x,beta))}
    iG = timesMatrix(uu,1/(lambda*delta+1),uu,n,n)
    U = iG%*%re
    if(length(m)>0){
      C=k%*%iG%*%U
      hat = rep(0,N)
      hat[m] = C
      hat[-m] = U
      U = hat
      }
    U = crossprod(U,Z)
    REML = list("VC"=VC,"Fixed"=B,"EBV"=U)
    
  }else{
    
    N = nrow(y)
    # Dealing with fixed effect matrix
    if(is.null(X)){X = matrix(1,N,1)}else{
      if(is.matrix(X)){if(nrow(X)!=N) stop("Fixed effect does not match dimensions of response variable")
      }else{if(class(X)=="formula"){X=model.matrix(X)}}}
    # Dealing with random effect
    if(is.null(K)&is.null(Z))stop("Either Z or K must be specified")
    if(is.null(K)){
      if(class(Z)=="formula"){Z=model.matrix(Z)-1}
      V=tcrossprod(Z)}
    if(is.null(Z)){V=K}
    if(is.null(Z)!=T&&is.null(K)!=T){
      if(class(Z)=="formula"){Z=model.matrix(Z)-1}
      V=crossprod(t(Z),K);V=tcrossprod(V,Z)};
    Z=diag(N); K=V
    
   ECM=function(Y,X,Z,K){Y=t(Y);X=t(X)
   ECM1=function(ytl, xtl, Vgt,Vet,Bt, deltal){Vlt=deltal*Vgt+Vet; invVlt=solve(Vlt+diag(1e-06,d))
   return(list(Vlt=Vlt, gtl=deltal*Vgt%*%invVlt%*%(ytl-Bt%*%xtl), Sigmalt=deltal*Vgt-deltal*Vgt%*%invVlt%*%(deltal*Vgt)))}
   wrapperECM1=function(l){ytl=Yt[,l]; xtl=Xt[,l]; deltal=eigZKZt$values[l]
   return( ECM1(ytl=ytl, xtl=xtl, Vgt=Vgt,Vet=Vet,Bt=Bt, deltal=deltal))}
   Vgfunc=function(l){Vgl=tcrossprod(outfromECM1[[l]]$gtl)
   return((1/n)*(1/eigZKZt$values[l])*(Vgl + outfromECM1[[l]]$Sigmalt))}
   Vefunc=function(l){etl = Yt[,l] - Bt%*%Xt[,l] - outfromECM1[[l]]$gtl
   return((1/n)*((tcrossprod(etl)+ outfromECM1[[l]]$Sigmalt)))}
   if (sum(is.na(Y))==0){ N=nrow(K); KZt=tcrossprod(K,Z)
   ZKZt=Z%*%KZt; eigZKZt = eigen(ZKZt); n=nrow(ZKZt); d=nrow(Y)
   Yt = Y%*%eigZKZt$vectors; Xt = X%*%eigZKZt$vectors
   Vgt =cov(t(Y))/2; Vet =cov(t(Y))/2; XttinvXtXtt=t(Xt)%*%solve(tcrossprod(Xt))
   Bt=Yt%*%XttinvXtXtt; Vetm1=Vet
   repeat{ outfromECM1=lapply(1:n, wrapperECM1)
         Vetm1=Vet; Gt=sapply(outfromECM1, function(x) {cbind(x$gtl)})
         Bt = (Yt - Gt) %*% XttinvXtXtt
         listVgts = lapply(1:n,Vgfunc); Vgt=Reduce('+', listVgts)
         listVets = lapply(1:n,Vefunc); Vet=Reduce('+', listVets)
         convnum=abs(sum(diag(Vet - Vetm1)))/abs(sum(diag(Vetm1)))
         convcond=tryCatch({convnum<1e-06}, error=function(e){return(FALSE)})
         if(convcond){break}}
 HobsInv=solve(kronecker(ZKZt,Vgt)+kronecker(diag(n),Vet)+diag(1e-06,d*n))
 ehat=matrix(Y - Bt%*%X,ncol=1, byrow=F)
 HobsInve=HobsInv%*%ehat; varvecG=kronecker(K,Vgt)
 gpred=varvecG%*%(kronecker(t(Z),diag(d)))%*%HobsInve
 Gpred=matrix(gpred, nrow=nrow(Y), byrow=F); 
 Fx=t(Bt); #rownames(Fx)=paste("beta",0:(nrow(Bt)-1),sep="")
 EBV = t(Gpred); colnames(Fx)=colnames(EBV)=paste("trait",1:d,sep="")
 VC = list("Vg"=Vgt,"Ve"=Vet,"h2"=diag(Vgt/(Vgt+Vet)),"GenCor"=cov2cor(Vgt))
 return(list("Fixed"=Fx,"VC"=VC,"EBV"=EBV))}}
 REML = ECM(Y=y,X=X,Z=Z,K=K)}

 class(REML) = "reml"
return(REML)}