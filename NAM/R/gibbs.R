
gibbs = function(y,Z=NULL,X=NULL,iK=NULL,iR=NULL,Iter=1500,Burn=500,Thin=4,DF=5,S=1,GSRU=FALSE){
  anyNA = function(x) any(is.na(x))  
  # Default for X; changing X to matrix if it is a formulas
  VY = var(y,na.rm=T)
  if(is.null(X)) X=matrix(1,length(y),1)
  if(class(X)=="formula"){
    X=model.frame(X)
    Fixes=ncol(X)
    XX=matrix(1,length(y),1)
    for(var in 1:Fixes) XX=cbind(XX,model.matrix(~X[,var]-1))
    X=XX
    rm(XX,Fixes)
  }
  
  # Defaults of Z: making "NULL","formula" and "matrix" as "list"
  if(is.null(Z)&is.null(iK)) stop("Either Z or iK must be specified")
  if(is.null(Z)) Z=list(diag(length(y)))
  if(class(Z)=="matrix") Z = list(Z)
  if(class(Z)=="formula") {
    Z=model.frame(Z)
    Randoms=ncol(Z)
    ZZ=list()
    for(var in 1:Randoms) ZZ[[var]]=model.matrix(~Z[,var]-1)
    Z=ZZ
    rm(ZZ,Randoms)
  }

  if(!GSRU){
    # Defaults for null and incomplete iK
    if(is.null(iK)){
      iK=list()
      Randoms=length(Z)
      for(var in 1:Randoms) iK[[var]]=diag(ncol(Z[[var]]))
    }
    if(class(iK)=="matrix") iK=list(iK)
    if(length(Z)!=length(iK)){
      a=length(Z)
      b=length(iK)
      if(a>b) for(K in 1:(a-b)) iK[[(K+b)]]=diag(ncol(Z[[K]]))
      if(b>a) for(K in 1:(b-a)) Z[[(K+a)]]=diag(ncol(iK[[K]]))
      rm(a,b)
    } 
  }
  
  # Predictiors should not have missing values
  if(any(is.na(X))|any(is.na(unlist(Z)))) stop("Predictors with missing values not allowed")
  
  # Add logit link function
  # Add residual variation
  
  # Gibb Sampling 
  # Garcia-Cortes, L. A. & Sorensen, D. (1996). GSE, 28(1), 121-126.
  # Sorensen, D., & Gianola, D. (2002). Springer.
  # Prior solution: de los Campos et al (2013). Genetics, 193(2), 327-345.
  
  # Thinning - which Markov Chains are going to be stored
  THIN = seq(Burn,Iter,Thin)
  
  # Some parameters with notation from the book
  nx = ncol(X)
  Randoms = length(Z) # number of random variables
  q = rep(0,Randoms); for(i in 1:Randoms) q[i]=ncol(Z[[i]])  
  N = nx+sum(q)
  
  # Qs1 and Qs2 regard where each random variable starts and ends, respectively
  Qs0 = c(nx,q)
  Variables = length(Qs0)
  Qs1 = Qs2 = rep(0,Variables)
  for(i in 1:Variables){
    Qs1[i]=max(Qs2)+1
    Qs2[i]=Qs1[i]+Qs0[i]-1
  }
  # Starting values for the variance components
  Ve = 1
  Va = lambda = rep(1,Randoms)
  
  # Linear system described as: WW+Sigma = Cg = r
  W = X
  for(i in 1:Randoms) W=cbind(W,Z[[i]])
  # MISSING
  W1=W
  if(any(is.na(y))){
    MIS = which(is.na(y))
    W=W[-MIS,]
    y=y[-MIS]
    if(!is.null(iR)) iR=iR[-MIS,-MIS]
  }
  n = length(y)
  # Keeping on
  # GSRU does not deal with WW or iK
  if(!GSRU){
  if(is.null(iR)){
    r = crossprod(W,y)
    WW = (crossprod(W))
  }else{
    r = crossprod(W,iR)%*%y
    WW = (crossprod(W,iR))%*%W
  }
  # Covariance Matrix
  Sigma = matrix(0,N,N)
  for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
  # Matching WW and Sigma
  C = WW+Sigma
  }else{
    L = rep(0,ncol(W))
    for(i in 1:Randoms) L[Qs1[i+1]:Qs2[i+1]] = lambda[i]
    xx = colSums(W^2)
  }
  g = rep(0,N)
  # Saving space for the posterior
  include = 0
  POSTg = matrix(0,ncol = length(THIN), N)
  POSTv = matrix(0,ncol = length(THIN), nrow = (Randoms+1))
  
  # Hperpriors: Degrees of freedom (DF) and Shape (S)
  df0 = DF
  
  if(is.null(S)){
    S0=rep(0,Randoms) # S0 has analystical solution (De los Campos et al 2013)
    if(!GSRU){
    for(i in 1:Randoms) 
      S0[i]=(var(y,na.rm=T)*0.5)/mean(
        (t((t(C[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]])-
              colMeans(C[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]]))))^2 )
    }else{
      for(i in 1:Randoms) 
        S0[i]=(var(y,na.rm=T)*0.5)/
        mean( (t((t(W[,Qs1[i+1]:Qs2[i+1]])-colMeans(W[,Qs1[i+1]:Qs2[i+1]]))))^2 )
    }
    
  }else{S0=rep(S,Randoms)}
    
  # Saving memory for some vectors
  e = rep(0,N)
  
  # Progression Bar
  pb=txtProgressBar(style=3)
  
  # LOOP
  for(iteration in 1:Iter){
    
    # Sampling hyper-priors
    S0a = runif(Randoms,S0*0.5,S0*1.5)
    df0a = runif(Randoms,df0*0.5,df0*1.5)
    dfu = q + df0a
    if(is.null(S)){ S0b = runif(1,0.0001,5) }else{ S0b = runif(1,S*0.5,S*1.5) }
    df0b = runif(1,min(2,df0)*0.5,min(2,df0)*1.5)
    dfe = n + df0b
    
    # Random variance
    for(i in 1:Randoms){
      # (ZiAZ+S0v0)/x2(v)
      if(!GSRU){
      Va[i] = (sum(crossprod(g[Qs1[i+1]:Qs2[i+1]],iK[[i]])*(g[Qs1[i+1]:Qs2[i+1]]))+  
                 S0a[i]*df0a[i])/rchisq(1,df=dfu[i])
      }else{
        Va[i] = (sum(crossprod(g[Qs1[i+1]:Qs2[i+1]]))+S0a[i]*df0a[i])/rchisq(1,df=dfu[i])
      }
    }
      
    # Residual variance
    e = y - W%*%g
    Ve = (crossprod(e)+S0b*df0b) / rchisq(1,df=dfe)
    
    # Ve/Va
    lambda = Ve/Va
    
    # Updating C
    if(!GSRU){
    for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
    C = WW+Sigma
    }
    
    # the C++ SAMP updates "g" and doesn't return anything
    if(!GSRU){
      SAMP(C,g,r,N,Ve)
    } else {
      SAMP2(W,g,y,xx,e,L,N,Ve) 
    }
        
    if(is.element(iteration,THIN)){
      include = include + 1;
      POSTg[,include] = g;
      POSTv[1:Randoms,include] = Va;
      POSTv[(Randoms+1),include] = Ve;
    }  
    # Advance progression bar
    setTxtProgressBar(pb,iteration/Iter)
  }
  # End progression bar
  close(pb)
  
  # Mode Function - Venter J.H. (1967). Ann. Math. Statist., 38(5):1446-1455.
  moda=function (x){
    it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
      y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
      diffs=sup-inf; i=min(which(diffs==min(diffs)))
      M=median(y[i:(i+k)]); it=it-1}; return(M)}
  
  rownames(POSTg)=paste("b",0:(N-1),sep="")
  for(i in 1:Randoms) rownames(POSTg)[Qs1[i+1]:Qs2[i+1]] = paste("u",i,".",1:Qs0[i+1],sep="")
  
  # Mean and Mode Posterior
  Mean.B = apply(POSTg,1,mean)
  Post.VC = c(apply(POSTv,1,moda))
  names(Post.VC) = c(paste("Va",1:Randoms,sep=""),"Ve")
  rownames(POSTv) = c(paste("Va",1:Randoms,sep=""),"Ve")
  
  # List of Coefficients
  Coefficients = list()
  Coefficients[[1]] = POSTg[1:nx,]
  for(i in 1:Randoms) Coefficients[[i+1]]=POSTg[(Qs1[i+1]:Qs2[i+1]),]
  names(Coefficients)[1]="Fixed"
  for(i in 1:Randoms) names(Coefficients)[i+1]=paste("Random",i,sep="")
  
  
  RESULTS = list(
                 "Coef.estimate" = Mean.B,
                 "VC.estimate" = Post.VC,
                 "Posterior.Coef" = Coefficients,
                 "Posterior.VC" = POSTv,
                 "Fit.mean" = W1%*%Mean.B
                 )
  
  class(RESULTS) = "gibbs"
  
  # Return
  return( RESULTS )
  
}

ml = function(y,Z=NULL,X=NULL,iK=NULL,iR=NULL,DF=-2,S=0){
  
  anyNA = function(x) any(is.na(x))
  
  # Default for X; changing X to matrix if it is a formulas
  VY = var(y,na.rm=T)
  if(is.null(X)) X=matrix(1,length(y),1)
  if(class(X)=="formula"){
    X=model.frame(X)
    Fixes=ncol(X)
    XX=matrix(1,length(y),1)
    for(var in 1:Fixes) XX=cbind(XX,model.matrix(~X[,var]-1))
    X=XX
    rm(XX,Fixes)
  }
  
  # Defaults of Z: making "NULL","formula" and "matrix" as "list"
  if(is.null(Z)&is.null(iK)) stop("Either Z or iK must be specified")
  if(is.null(Z)) Z=list(diag(length(y)))
  if(class(Z)=="matrix") Z = list(Z)
  if(class(Z)=="formula"){
    Z=model.frame(Z)
    Randoms=ncol(Z)
    ZZ=list()
    for(var in 1:Randoms) ZZ[[var]]=model.matrix(~Z[,var]-1)
    Z=ZZ
    rm(ZZ,Randoms)
  }
  
  # Defaults for null and incomplete iK
  if(is.null(iK)){
    iK=list()
    Randoms=length(Z)
    for(var in 1:Randoms) iK[[var]]=diag(ncol(Z[[var]]))
  }
  if(class(iK)=="matrix") iK=list(iK)
  if(length(Z)!=length(iK)){
    a=length(Z)
    b=length(iK)
    if(a>b) for(K in 1:(a-b)) iK[[(K+b)]]=diag(ncol(Z[[K]]))
    if(b>a) for(K in 1:(b-a)) Z[[(K+a)]]=diag(ncol(iK[[K]]))
    rm(a,b)
  }
  
  # Predictiors should not have missing values
  if(any(is.na(X))|any(is.na(unlist(Z)))) stop("Predictors with missing values not allowed")
  
  # Some parameters with notation from the book
  nx = ncol(X)
  Randoms = length(Z) # number of random variables
  q = rep(0,Randoms); for(i in 1:Randoms) q[i]=ncol(Z[[i]])  
  N = nx+sum(q)
  
  # Qs1 and Qs2 regard where each random variable starts and ends, respectively
  Qs0 = c(nx,q)
  Variables = length(Qs0)
  Qs1 = Qs2 = rep(0,Variables)
  for(i in 1:Variables){
    Qs1[i]=max(Qs2)+1
    Qs2[i]=Qs1[i]+Qs0[i]-1
  }
  
  # Starting values for the variance components
  Ve = 1
  Va = lambda = rep(1,Randoms)
  
  # Linear system described as: WW+Sigma = Cg = r
  W = X
  for(i in 1:Randoms) W=cbind(W,Z[[i]])
  
  # MISSING
  W1=W
  if(any(is.na(y))){
    MIS = which(is.na(y))
    W=W[-MIS,]
    y=y[-MIS]
    if(!is.null(iR)) iR=iR[-MIS,-MIS]
  }
  n = length(y)
  
  if(is.null(iR)){
    r = crossprod(W,y)
    WW = (crossprod(W))
  }else{
    r = crossprod(W,iR)%*%y
    WW = (crossprod(W,iR))%*%W
  }
  
  # Covariance Matrix
  Sigma = matrix(0,N,N)
  for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
  # Matching WW and Sigma
  C = WW+Sigma
  
  g = rep(0,N)
  
  # Variance components
  dfu = q+DF
  dfe = n+DF
  
  # Saving memory for some vectors
  e = rep(0,N)
  
  # convergence factor
  VC = c(Va,Ve)
  cf = 1
  
  # LOOP
  while(cf>1e-8){
    
    # Ve/Va
    lambda = Ve/Va
    
    # Updating C
    for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = iK[[i]]*lambda[i]
    C = WW+Sigma
    
    # Residual variance
    e = y - tcrossprod(g,W)
    Ve = (tcrossprod(e)+S*DF)/dfe
    
    # the C++ SAMP updates "g" and doesn't return anything
    gs(C,g,r,N)
    
    # Random variance
    for(i in 1:Randoms){
      SS = S*DF+(sum(crossprod(g[Qs1[i+1]:Qs2[i+1]],iK[[i]])*(g[Qs1[i+1]:Qs2[i+1]])))
      SS = max(SS,1e-8)
      Va[i] = SS/dfu[i]
    }
    
    # convergence
    cf = sum((c(Va,Ve)-VC)^2)
    VC = c(Va,Ve)
    
  }
  
  names(g) = paste("b",0:(N-1),sep="")
  for(i in 1:Randoms) names(g)[Qs1[i+1]:Qs2[i+1]] = paste("u",i,".",1:Qs0[i+1],sep="")
  names(VC) = c(paste("Va",1:Randoms,sep=""),"Ve")
  
  # List of Coefficients
  
  RESULTS = list(
    "Coef" = g,
    "VC" = VC,
    "Fit" = W1%*%g
  )
  
  # Return
  return( RESULTS )
  
}

plot.gibbs = function(x,...){
  anyNA = function(x) any(is.na(x))
  par(ask=TRUE)
  vc = nrow(x$Posterior.VC)-1
  for(i in 1:vc) plot(density(x$Posterior.VC[i,],...),main=paste("Posterior: Term",i,"variance"))
  plot(density(x$Posterior.VC[vc+1,]),main=paste("Posterior: Residual Variance"),...)
  par(ask=FALSE)
}


gibbs2 = function(Y,Z=NULL,X=NULL,iK=NULL,Iter=150,Burn=50,Thin=3,DF=5,S=1){
  
  anyNA = function(x) any(is.na(x))  
  
  Y0 = Y
  Q = ncol(Y)
  n0 = nrow(Y)
  mNa = !is.na(Y)
  m0 = crossprod(mNa)
  eAdj = n0/m0
  E = matrix(0,n0,Q)
  
  # Default for X; changing X to matrix if it is a formulas
  VY = apply(Y,2,var,na.rm=T)
  if(is.null(X)) X=matrix(1,nrow(Y),1)
  if(class(X)=="formula"){
    X=model.frame(X)
    Fixes=ncol(X)
    XX=matrix(1,n0,1)
    for(var in 1:Fixes) XX=cbind(XX,model.matrix(~X[,var]-1))
    X=XX
    rm(XX,Fixes)
  }
  
  # Defaults of Z: making "NULL","formula" and "matrix" as "list"
  if(is.null(Z)&is.null(iK)) stop("Either Z or iK must be specified")
  if(is.null(Z)) Z=list(diag(n0))
  if(class(Z)=="matrix") Z = list(Z)
  if(class(Z)=="formula") {
    Z=model.frame(Z)
    Randoms=ncol(Z)
    ZZ=list()
    for(var in 1:Randoms) ZZ[[var]]=model.matrix(~Z[,var]-1)
    Z=ZZ
    rm(ZZ,Randoms)
  }
  
  if(is.null(iK)){
    iK=list()
    Randoms=length(Z)
    for(var in 1:Randoms) iK[[var]]=diag(ncol(Z[[var]]))
  }  
  if(class(iK)=="matrix") iK=list(iK)
  if(length(Z)!=length(iK)){
    a=length(Z)
    b=length(iK)
    if(a>b) for(K in 1:(a-b)) iK[[(K+b)]]=diag(ncol(Z[[K]]))
    if(b>a) for(K in 1:(b-a)) Z[[(K+a)]]=diag(ncol(iK[[K]]))
    rm(a,b)
  }
  
  # Predictiors should not have missing values
  if(any(is.na(X))|any(is.na(unlist(Z)))) stop("Predictors with missing values not allowed")
  
  # Thinning - which Markov Chains are going to be stored
  THIN = seq(Burn,Iter,Thin)
  
  # Some parameters with notation from the book
  # Adapted for Multi-trait
  nx = ncol(X)*Q
  Randoms = length(Z) # number of random variables
  MSx = rep(0,Randoms)
  q = rep(0,Randoms); for(i in 1:Randoms){
    q[i]=ncol(Z[[i]])
    MSx[i]=mean(colSums(Z[[i]]^2))
  }
  q = q*Q
  N = nx+sum(q)
  
  # Qs1 and Qs2 regard where each random variable starts and ends, respectively
  Qs0 = c(nx,q)
  Variables = length(Qs0)
  Qs1 = Qs2 = rep(0,Variables)
  for(i in 1:Variables){
    Qs1[i]=max(Qs2)+1
    Qs2[i]=Qs1[i]+Qs0[i]-1
  }
  
  # Priors
  Sp = ifelse(!is.null(S),S,0.5*VY*(DF+2)/MSx)
  S0 = diag(Sp*DF,Q)
  df0a = DF+q/Q
  df0e = DF+n0  
  
  # Starting values for the variance components
  Ve = 0.1*diag(diag(var(Y,na.rm = T)))+1e-4
  Va = lambda = list()
  for(i in 1:Randoms){
    Va[[i]] = 0.01+diag(VY)+1e-4
    lambda[[i]] = solve(Va[[i]])
  }
  
  # KRONECKERS
  Yk=matrix(Y)
  Wk = kronecker(diag(Q),X)
  for(i in 1:Randoms) Wk=cbind(Wk,kronecker(diag(Q),Z[[i]]))
  R = kronecker(Ve,diag(n0))
  
  # MISSING
  W1=Wk
  if(any(is.na(Yk))){
    Ms = TRUE
    MIS = which(is.na(Yk))
    Wk=Wk[-MIS,]
    Yk=Yk[-MIS]
    R = R[-MIS,-MIS]
  }else{
    Ms = FALSE
  }
  n = length(Yk)
  
  # Keeping on
  iR = chol2inv(R)
  MM = t(Wk) %*% iR %*% Wk
  r = t(Wk) %*% iR %*% Yk
  
  # Covariance Matrix
  Sigma = matrix(0,N,N)
  for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = kronecker(lambda[[i]],iK[[i]])
  
  # Matching WW and Sigma
  C = MM+Sigma
  g = rep(0,N)
  
  # Saving space for the posterior
  include = 0
  POSTa = array(data = 0,dim = c(Q,Q,Randoms,length(THIN)))
  POSTe = array(data = 0,dim = c(Q,Q,length(THIN)))
  POSTg = matrix(0,N,length(THIN))
  
  # Saving memory for some vectors
  e = rep(0,n)
  
  # Progression Bar
  pb=txtProgressBar(style=3)
  
  # LOOP
  for(iteration in 1:Iter){
    
    # the C++ SAMP updates "g" and doesn't return anything
    SAMP(C,g,r,N,1)
    
    # Residual variance
    e[1:n] = Yk - Wk%*%g
    E[mNa] = e
    SSe = eAdj*crossprod(E) + S0
    #E2 = tapply(e,obs_index,crossprod)
    Ve = solve(rWishart(1,df0e,solve(SSe))[,,1])
    
    # Random variance
    for(i in 1:Randoms){
      u = matrix(g[Qs1[i+1]:Qs2[i+1]],ncol = Q)
      SSa = S0+t(u)%*%iK[[i]]%*%u
      lambda[[i]] = rWishart(1,df0a[i],solve(SSa))[,,1]
      Va[[i]] = solve(lambda[[i]])
    }
    
    # Updating C
    
    R = kronecker(Ve,diag(n0))
    if(Ms) R=R[-MIS,-MIS]
    
    iR = chol2inv(R)
    MM = t(Wk) %*% iR %*% Wk
    r = t(Wk) %*% iR %*% Yk
    for(i in 1:Randoms) Sigma[Qs1[i+1]:Qs2[i+1],Qs1[i+1]:Qs2[i+1]] = kronecker(lambda[[i]],iK[[i]])
    C = MM+Sigma
    
    # Storing elements into posteriors
    if(is.element(iteration,THIN)){
      include = include + 1
      POSTg[,include] = g
      # Random: Trait,Trait,Random variable,Iteration
      for(i in 1:Randoms) POSTa[,,i,include] = Va[[i]]
      # Residual: Trait,Trait,Iteration
      POSTe[,,include] = Ve
    }  
    # Advance progression bar
    setTxtProgressBar(pb,iteration/Iter)
  }
  # End progression bar
  close(pb)
  
  rownames(POSTg)=1:N
  namesX = paste(rep(paste('b',0:(nx/2-1),sep=''),Q),sort(rep(1:Q,nx/2)),sep='_trait')
  rownames(POSTg)[1:(nx)] = namesX
  
  for(i in 1:Randoms){
    LZ = length(Qs1[i+1]:Qs2[i+1]) # length of Z_i
    NZ1 = paste(paste('u',1:(LZ/2),sep=''),'_eff',i,sep='')
    NZ2 = paste('trait',sort(rep(1:Q,LZ/2)),sep='')
    NZ3 = rep(NZ1,length.out=length(NZ2))
    namesZ = paste(NZ3,NZ2,sep='_')
    rownames(POSTg)[Qs1[i+1]:Qs2[i+1]] = namesZ
  } 
  
  #####################
  ###               ###
  ### NEW PROBLEMS  ###
  ###               ###
  #####################
  
  # Mean and Mode Posterior
  Coef = rowMeans(POSTg)
  
  VCE = Ve
  for(i in 1:Q){
    for(j in 1:Q){
      VCE[i,j] = mean(POSTe[i,j,])
    }}
  
  VCA = Va
  for(i in 1:Q){
    for(j in 1:Q){
      for(k in 1:Randoms){
        VCA[[k]][i,j] = mean(POSTa[i,j,k,])
      }}}
  
  
  namesVCs = paste('trait',1:Q,sep='')
  namesVCs = list(namesVCs,namesVCs)
  dimnames(VCE) = namesVCs
  names(VCA) = paste('term',1:Randoms,sep='')
  for(i in 1:Randoms) dimnames(VCA[[i]]) = namesVCs
  
  
  # Output
  Posterior = list('Coef'=POSTg,'VarA'=POSTa,'VarE'=POSTe)
  HAT = matrix(W1%*%Coef,ncol = Q)
  
  RESULTS = list(
    "Coef" = Coef,
    "VarA" = VCA,
    "VarE" = VCE,
    "Posterior" = Posterior,
    "Fit" = HAT
  )
  
  # Return
  return( RESULTS )
  
}

covar = function(sp=NULL,rho=3.5,type=1,dist=2.5){
  if(is.null(sp)) {
    sp = cbind(rep(1,49),rep(c(1:7),7),as.vector(matrix(rep(c(1:7),7),7,7,byrow=T)))
    colnames(sp)=c("block","row","col")
    cat("Example of field information input 'sp'\n")
    print(head(sp,10))
  }
  if(type==1) cat("Exponential Kernel\n")
  if(type==2) cat("Gaussian Kernel\n")
  obs=nrow(sp)
  sp[,2]=sp[,2]*dist
  fields=sp[,1]
  Nfield=length(unique(fields))
  quad=matrix(0,obs,obs)
  for(j in 1:Nfield){
    f=which(fields==j);q=sp[f,2:3]    
    e=dist(q);e=as.matrix(e);e=round(e,3)
    d=exp(-e^type/(rho^type))
    d2=round(d,2);quad[f,f]=d2} 
  if(obs==49){
    M=matrix(quad[,25],7,7)
    dimnames(M) = list(abs(-3:3),abs(-3:3))
    print(M)
  } 
  if(obs!=49) return(quad)
}

PedMat = function(ped=NULL){
  if(is.null(ped)){
    id = 1:11
    dam = c(0,1,1,1,1,2,0,4,6,8,9)
    sire = c(0,0,0,0,0,3,3,5,7,3,10)
    example = cbind(id,dam,sire)
    cat('Example of pedigree\n')
    print(example)
    cat('- It must follow chronological order\n')
    cat('- Zeros are used for unknown\n')
  }else{
    n = nrow(ped)
    A=diag(0,n)
    for(i in 1:n){
      for(j in 1:n){
        if(i>j){A[i,j]=A[j,i]}else{
          d = ped[j,2]
          s = ped[j,3]
          if(d==0){Aid=0}else{Aid=A[i,d]}
          if(s==0){Ais=0}else{Ais=A[i,s]}
          if(d==0|s==0){Asd=0}else{Asd=A[d,s]}
          if(i==j){Aij=1+0.5*Asd}
          if(i!=j){Aij=0.5*(Aid+Ais)}
          A[i,j]=Aij}}}
    return(A)}}

PedMat2 = function (ped,gen=NULL,IgnoreInbr=FALSE,PureLines=FALSE){
  
  n = nrow(ped)
  A = diag(0, n)
  
  if(is.null(gen)){
    
    # WITHOUT GENOTYPES
    for (i in 1:n) {
      for (j in 1:n) {
        if (i > j) {
          A[i, j] = A[j, i]
        } else {
          d = ped[j, 2]
          s = ped[j, 3]
          if (d == 0) Aid = 0 else Aid = A[i, d]
          if (s == 0) Ais = 0 else Ais = A[i, s]
          if (d == 0 | s == 0) Asd = 0 else Asd = A[d, s]
          if (i == j) Aij = 1 + 0.5 * Asd else Aij = 0.5 * (Aid + Ais)
          A[i, j] = Aij
        }}}
    
  }else{
    
    # WITH GENOTYPES 
    G = as.numeric(rownames(gen))
    if(any(is.na(gen))){
      ND = function(g1,g2){
        X = abs(g1-g2)
        L = 2*sum(!is.na(X))
        X = sum(X,na.rm=TRUE)/L
        return(X)} 
    }else{
      ND = function(g1,g2) sum(abs(g1-g2))/(2*length(g1))
    }
    Inbr = function(g) 2-mean(g==1)
    
    # LOOP     
    for (i in 1:n) {
      for (j in 1:n) {
        
        #######################
        if (i > j) {
          A[i, j] = A[j, i]
        } else {
          
          d = ped[j, 2]
          s = ped[j, 3]
          
          if(j%in%G){
            
            if(d!=0&d%in%G){ Aid=2*ND(gen[paste(j),],gen[paste(d),]) }else{if(d==0) Aid=0 else Aid=A[i,d]}
            if(s!=0&s%in%G){ Ais=2*ND(gen[paste(j),],gen[paste(s),]) }else{if(s==0) Ais=0 else Ais=A[i,s]}
            Asd = Inbr(gen[paste(j),])
            
          }else{
            if (d == 0) Aid = 0 else Aid = A[i, d]
            if (s == 0) Ais = 0 else Ais = A[i, s]
            if (d == 0 | s == 0) Asd = 0 else Asd = A[d, s]
          }
          
          if (i == j) Aij = ifelse(IgnoreInbr,1,ifelse(PureLines,ifelse(i%in%G,Asd,2),1+0.5*Asd))
          
          else Aij = 0.5*(Aid+Ais)
          A[i, j] = round(Aij,6)
        }
        #######################
        
      }}
  }
  return(A)
}
