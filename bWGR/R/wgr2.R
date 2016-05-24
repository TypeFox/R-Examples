wgr2 = function(Y,gen,it=1000,bi=250,th=3,df=5,eigK=NULL,EigT=0.05,R2=0.5,pi=0,verb=FALSE){
  
  Poly = 0
  if(!is.null(eigK)) Poly = TRUE
  X = gen
  xx = colSums(X^2)
  MSx = sum(apply(X,2,var,na.rm=T))
  S_prior = R2*apply(Y,2,var,na.rm=T)*(df+2)/MSx/(1-pi)
  k = ncol(Y)
  n = nrow(Y)
  p = ncol(X)
  s = diag(S_prior,k)
  N = crossprod(!is.na(Y))
  mN = mean(N)
  B = G = D = matrix(0,p,k)
  E = apply(Y,2,function(x)x-mean(x,na.rm = T))
  E[is.na(E)]=0
  md = rep(pi,k)
  mu = colMeans(Y,na.rm = T)
  v = var(Y,na.rm = T)
  VA = 0.4*diag(diag(var(Y,na.rm = T)))
  VE = 0.1*diag(diag(var(Y,na.rm = T)))
  
  # Polygenic term
  if(Poly){
    U = eigK$vectors
    V = eigK$values
    pk = sum(V>EigT)
    if(!is.null(EigT)) {U=U[,1:pk];V=V[1:pk]}
    sk_prior = R2*apply(Y,2,var,na.rm=T)*(df+2)/pk
    s2 = diag(sk_prior*df,k)
    xx2 = rep(1,p)
    VK = 0.1*diag(diag(var(Y,na.rm = T)))
    mcP = P = matrix(0,pk,k)
    mcVK = matrix(0,k,k)
  }  
  
  # Store
  mc = seq(bi,it,th)
  lmc = length(mc)
  mcMu = rep(0,k)
  mcB = mcG = mcD = B
  mcVA = mcVE = matrix(0,k,k)
  
  # Indicators
  y = z = list()
  for(i in 1:k){
    z[[i]] = which(!is.na(Y[,i]))
    y[[i]] = Y[z[[i]],i]
  }
  
  # MCMC
  if(verb) pb = txtProgressBar(style=3)
  
  for(j in 1:it){
    
    mEA = diag(solve(VA)%*%VE)
    if(Poly) mEK = diag(solve(VK)%*%VE)
    if(pi>0){PI = rbeta(1,10*pi+md+1,10*(1-pi)-md+1)}else{PI=0}
    
    # Update regression coefficients
    for(i in 1:k){
      
      lambda = mEA[i]
      L = rep(lambda,p)
      up = KMUP(X[z[[i]],],B[,i],xx,E[z[[i]],i],L,p,VE[i,i],PI)
      b = up[[1]]
      d = up[[2]]
      e = up[[3]]
      if(pi>0) d[is.nan(d)] = 1
      B[,i] = b
      D[,i] = d
      G[,i] = up[[1]]*d
      
      if(Poly){
        lambda2 = mEK[i]
        L2 = lambda2/V
        up = KMUP(U[z[[i]],],P[,i],xx2,e,L2,pk,VE[i,i],0)
        u = up[[1]]
        e = up[[3]]
        P[,i] = u
      } 
      
      mu[i] = rnorm(1,mu[i]+mean(e),VE[i,i]/N[i,i])
      E[z[[i]],i] = y[[i]]-X[z[[i]],]%*%B[,i]-mu[i]
      md[i] = mean(D[,i])
    }
    
    # Update variance components
    VA = solve(rWishart(1,n+df,solve(crossprod(B)+s))[,,1])
    if(Poly) VK = solve(rWishart(1,pk+df,solve(crossprod(P/V,P)+s2))[,,1])
    SSe = (crossprod(E)/N)*mN
    VE = solve(rWishart(1,mN+2,solve(SSe))[,,1])

    if(j%in%mc){
      mcMu = mcMu + mu
      mcB = mcB + B
      mcG = mcG + G
      mcD = mcD + D
      mcVA = mcVA + VA
      mcVE = mcVE + VE
      if(Poly){
        mcVK = mcVK + VK
        mcP = mcP + P
      } 
    }
    
    if(verb) setTxtProgressBar(pb, j/it)
  }
  
  if(verb) close(pb)
  
  # Posterior means
  mcMu = mcMu/lmc
  mcB = mcB/lmc
  mcG = mcG/lmc
  mcD = mcD/lmc
  mcVA = mcVA/lmc
  mcVE = mcVE/lmc
  
  if(Poly){
    mcP = mcP/lmc
    mcVK = mcVK/lmc
    polygenic = U%*%mcP
    A = t(mcMu+t(X%*%mcB))+polygenic
    final = list('Mu'=mcMu,'B'=mcB,'G'=mcG,'D'=mcD,'U'=polygenic,
                 'VA'=mcVA,'VK'=mcVK,'VE'=mcVE,'Fit'=A)
  }else{
    A = t(mcMu+t(X%*%mcB))
    final = list('Mu'=mcMu,'B'=mcB,'G'=mcG,'D'=mcD,
                 'VA'=mcVA,'VE'=mcVE,'Fit'=A)
  }
  
  return(final)
}