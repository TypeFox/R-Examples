dcovustatC <- function( x, y, alpha=1)
  .Call("dcovustatC", x, y, alpha, PACKAGE = "steadyICA")


gradmdcov <- function (Z1, Z2) 
  .Call("gradmdcov", Z1, Z2, PACKAGE = "steadyICA")

givens.rotation <- function(theta=0, d=2, which=c(1,2))
{
  # David S. Matteson
  # 2008.04.28
  # For a given angle theta, returns a d x d Givens rotation matrix
  #
  # Ex: for i < j , d = 2:  (c -s)
  #                         (s  c)
  c = cos(theta)
  s = sin(theta)
  M = diag(d)
  a = which[1]
  b = which[2]
  M[a,a] =  c
  M[b,b] =  c
  M[a,b] = -s
  M[b,a] =  s
  M
}

theta2W = function(theta)
{
  # David S. Matteson
  # 2011.06.27
  # For a vector of angles theta, returns W, a d x d Givens rotation matrix:
  # W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2 
  ##  if(theta < 0  || pi < theta){stop("theta must be in the interval [0,pi]")}
  d = (sqrt(8*length(theta)+1)+1)/2
  if(d - floor(d) != 0){stop("theta must have length: d(d-1)/2")}
  W = diag(d)
  index = 1
  for(j in 1:(d-1)){
    for(i in (j+1):d){
      Q.ij = givens.rotation(theta[index], d, c(i,j))
      W = Q.ij %*% W 
      index = index + 1
    }
  }
  W
}

W2theta = function(W)
{
  # David S. Matteson
  # 2011.06.27
  # Decompose a d by d orthogonal matrix W into the product of
  # d(d-1)/2 Givens rotation matrices. Returns theta, the d(d-1)/2 by 1
  # vector of angles, theta.
  # W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2 
  if(dim(W)[1] != dim(W)[2]){stop("W must be a square matrix")}
  W = t(W) 
  d = dim(W)[1]
  #  if(sum(abs(t(W)%*%W  - diag(d))) > 1e-10){stop("W must be an orthogonal matrix")}
  theta = numeric(d*(d-1)/2)
  index = 1
  for(j in 1:(d-1)){
    for(i in (j+1):d){      
      x = W[j,j]
      y = W[i,j]        
      theta.temp = atan2(y,x)    
      Q.ij = givens.rotation(theta.temp, d, c(i,j))
      W.temp = Q.ij %*% W  
      W = W.temp
      theta[index] = theta.temp
      index = index + 1
    }
  }
  theta 
}

dcovICA <- function(Z, theta.0 = 0){
  fun.theta = function(theta){
    W = theta2W(theta)
    SH = (Z %*% t(W)) 
    dcovustat(SH[,1], SH[,2]) 
  }  
  out = optimize(f = fun.theta, interval = c(theta.0, theta.0 + pi/2))
  theta.hat = out$minimum
  W.hat=theta2W(theta.hat)
  #W.hat saved for parameterization X = S A
  #S saved for parameterization X = S A
  list(theta.hat = theta.hat, W = t(W.hat), S=Z%*%t(W.hat),obj = out$objective)
}

# Wrapper for c function
dcovustat = function(x,y,alpha=1){
  x = as.matrix(x)
  y = as.matrix(y)
  x = x/nrow(x)
  y = y/nrow(y)
  dcovustatC(x,y,alpha)
}

# Calculate either the symmetric and asymmetric 
# multivariate dcov statistic
multidcov <- function(S,symmetric=TRUE,alpha=1) {
  d = ncol(S)
  md = 0
  if(symmetric)  {
    for (i in 1:d)  md = md + dcovustat(S[,i],S[,-i],alpha)
  }    else {
    for (i in 1:(d-1)) md = md + dcovustat(S[,i],S[,(i+1):d],alpha)
  }
  md
}

# Asymmetric multivariate dcov
#adcov <- function(S) {
#  d = ncol(S)
#  md = 0
#  for (i in 1:(d-1)) {
#    md = md + dcovustat(S[,i],S[,(i+1):d])
#  }
#  md
#}

gradCpp = function(q,S) {
  Z1 = S[,q]/nrow(S)
  Z2 = as.matrix(S[,-q]/nrow(S))
  temp = gradmdcov(Z1,Z2)
  gS = matrix(0,nrow=nrow(S),ncol=ncol(S))
  gS[,q]=temp$gSq
  gS[,-q]=temp$gSnq
  gS
}


#-----------------------------------
# R version for systems without RCPP.
# This is MUCH slower:
gradR<-function(q,S) {
  #calculates gradient of dcov of S_q (n x 1 vector) and S_nq (n x (d-1))
  #note that first rv is univariate.
  d=ncol(S)
  n=nrow(S)
  #objects needed for both grad S_q and grad S_nq:
  S.nq=S[,-q]
  dist.S.nq=as.matrix(dist(S.nq))
  coef.T1=2/(n*(n-1))
  coef.T2=2/(n*(n-1))
  coef.T3=2/(n*(n-1)*(n-2))
  colSums.dist.S.nq=colSums(dist.S.nq*coef.T3) #*coef.T3 to prevent overflow
  
  #----------------
  #grad S_q:
  temp.S.q=matrix(S[,q],nrow=n,ncol=n)
  diff.S.q=t(temp.S.q)-temp.S.q #n x n anti-symmetric matrix. Columns 
  #correspond to S[i,q]-S[,q]
  rm(temp.S.q)
  g.dist.S.q=-1L*(diff.S.q<0)+1L*(diff.S.q>0) #Columns correspond to the derivative of the distance between scalars.
  dist.S.q=abs(diff.S.q) #Used in s_p terms
  rm(diff.S.q)
  
  #T1 S_q terms:
  g.T1.q=g.dist.S.q*dist.S.nq
  g.T1.q=colSums(g.T1.q*coef.T1) 
  
  #T2 S_q terms:
  T2.y=sum(colSums.dist.S.nq/(coef.T3*(n*(n-1))))  
  g.T2.x=apply(g.dist.S.q*coef.T2,2,sum)
  g.T2.q=g.T2.x*T2.y
  
  #T3 S_q terms:
  c.dist.s.nq=matrix(colSums.dist.S.nq,nrow=n,ncol=n)
  int.dist.s.nq=g.dist.S.q*(c.dist.s.nq+t(c.dist.s.nq))
  temp.g.T3.q=colSums(int.dist.s.nq)
  g.T3.q=temp.g.T3.q-2*g.T1.q*coef.T3/coef.T1
  rm(c.dist.s.nq,int.dist.s.nq,temp.g.T3.q)
  #grad S_q
  g.q=g.T1.q+g.T2.q-g.T3.q
  rm(g.T1.q,g.T2.q,g.T3.q)
  
  #------------------------
  # grad S_nq:
  r.dist.S.nq=1/dist.S.nq
  r.dist.S.nq[r.dist.S.nq==Inf]=0
  #dist.S.q=as.matrix(dist(S[,q]))
  colSums.dist.S.q=colSums(dist.S.q)
  c.dist.S.q=matrix(colSums.dist.S.q,nrow=n,ncol=n)
  int.dist.S.q=c.dist.S.q+t(c.dist.S.q)
  T2.x=sum(colSums.dist.S.q/(n*(n-1)))
  
  index.nq=c(1:d)[-q]
  g.T1.nq=matrix(0,n,d)
  g.T2.y=g.T1.nq
  g.T3.nq=g.T1.nq
  
  for(p in index.nq) {
    temp.S.p=matrix(S[,p],nrow=n,ncol=n)
    diff.S.p=t(temp.S.p)-temp.S.p #n x n anti-symmetric matrix.
    g.dist.S.p.nq=diff.S.p*r.dist.S.nq 
    g.T1.p=coef.T1*dist.S.q*g.dist.S.p.nq
    g.T1.nq[,p]=apply(g.T1.p,2,sum)
    g.T2.y[,p]=apply(coef.T2*g.dist.S.p.nq,2,sum)
    temp=g.dist.S.p.nq*int.dist.S.q
    temp.g.T3.p=colSums(coef.T3*temp)
    g.T3.nq[,p]=temp.g.T3.p-2*g.T1.nq[,p]*(coef.T3/coef.T1)    
  }
  #grad S_nq:
  g.nq=g.T1.nq+T2.x*g.T2.y-g.T3.nq
  
  #grad S:
  g.nq[,q]=g.q
  g.nq
}

#-------------------------------------------------------------
##For general use: Symmetric dCov
grad.mdCov=function(S,method = c('Cpp','R')) {
  method = match.arg(method)
  d=ncol(S)
  n=nrow(S)
  m.grad=matrix(0,n,d)
  if (method == 'Cpp') {
    for(i in 1:d) m.grad = m.grad + gradCpp(q=i,S=S)
  }
  else if (method == 'R') {
    for(i in 1:d) m.grad = m.grad + gradR(q=i,S=S)
  }
  m.grad
}
  
grad.adCov=function(S,method = c('Cpp','R')) {
   method = match.arg(method)
   d=ncol(S)
   n=nrow(S)
   m.grad=matrix(0,n,d)
   if (method == 'Cpp') {
     for(i in 1:(d-1)) m.grad[,i:d] = m.grad[,i:d] + gradCpp(q=1,S=S[,i:d])
   }
   else if (method == 'R') {
     for(i in 1:(d-1)) m.grad[,i:d] = m.grad[,i:d] + gradR(q=1,S=S[,i:d])
   }
  m.grad
}

#-------------------
# Whitening Function:
whitener <- function(X,n.comp=ncol(X),center.row=FALSE,irlba=FALSE) {
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Creates model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  if(irlba==FALSE) svd.x=svd(x.center,nu=n.comp,nv=n.comp)
  if(irlba==TRUE) {
    requireNamespace(irlba)
    svd.x=irlba(x.center,nu=n.comp,nv=n.comp)
  }
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)  
  return(list(whitener=t(ginv(svd.x$v%*%diag(svd.x$d[1:n.comp])/sqrt(n.rep-1))),Z=sqrt(n.rep-1)*svd.x$u,mean=apply(X,2,mean)))
}

# Benjamin Risk 8 November 2015: Changed name to steadyICA
steadyICA <- function(X, n.comp = ncol(X), w.init= NULL, PIT = FALSE, bw='SJ',adjust=1,whiten = FALSE, irlba = FALSE, symmetric=FALSE,eps = 1e-08, alpha.eps = 1e-08, maxit = 100, method = c('Cpp','R'), verbose = FALSE) {
  p <- ncol(X)
  d <- n.comp
  pd <- p * n.comp
  oldW=w.init
  if(p != d && PIT==FALSE && whiten==FALSE) stop("For n.comp != p, must use PIT or pre-whitening")
  if(is.null(oldW)){
    if (p == d | whiten==TRUE) {
      oldW = diag(n.comp)
  } else {
    oldW = matrix(c(diag(n.comp),rep(0,pd-d*d)),nrow=p,ncol=d,byrow=TRUE)
    }
  }
  if(whiten) {
    zData = whitener(X,n.comp=n.comp,irlba=irlba)
    whiteU <- zData$whitener
    Z = zData$Z
    rm(zData)
  } else {
    Z = X
  }
  if(p>nrow(Z)) warning("X must be n x p")
  if(n.comp==2 && PIT==FALSE) message('for n.comp=2, you should use multidcov::dcovICA -- it is much faster')  
  #if(sum(abs(cov(Z))) > n.comp+n.comp*0.005) warning("X must be pre-whitened -- check that covariance equals identity")
  method = match.arg(method)
  alpha <- 1
  S <- Z%*%oldW
  if(PIT) { 
    pS <- est.PIT(S=S,bw=bw,adjust=adjust)
    S <- pS$Fx*sqrt(12) #To have (asymptotically) unit variance; then alpha scales appropriately; 
    fS <- pS$fx*sqrt(12)
    rm(pS)
  } else {
    fS <- 1
  }
  
  if(symmetric) {
    deltaS <- grad.mdCov(S,method=method)*fS
  } else {
    deltaS <- grad.adCov(S,method=method)*fS
  }
  curF <- multidcov(S,symmetric=symmetric)
  deltaW <- crossprod(Z,deltaS)
  Table <- NULL
  iter <- 1
  deltaF <- -1
  while (iter < maxit) {
    alpha <- 2 * alpha
    tempW <- oldW - alpha * deltaW
    if(p==d | whiten==TRUE) {
      UDV <- svd(tempW)
      newW <- tcrossprod(UDV$u,UDV$v)
    } else {
      newW <- tempW
    }
    S <- Z %*% newW
    if(PIT) {
      pS <- est.PIT(S=S,bw=bw,adjust=adjust)
      S <- pS$Fx*sqrt(12) #To have (asymptotically) unit variance; then alpha scales appropriately; 
      fS <- pS$fx*sqrt(12)
      rm(pS)
    }
    
    newF <- multidcov(S=S,symmetric=symmetric) 
    while (newF >= curF) { #Search for alpha that reduces f
      alpha <- alpha/2
      if(alpha < alpha.eps) { 
        warning("alpha is less than alpha.eps -- if dcov is still changing, then try a different w.init")
        break
      }
      tempW <- oldW - alpha * deltaW
      if(p==d | whiten==TRUE) {
        UDV <- svd(tempW)
        newW <- tcrossprod(UDV$u,UDV$v)
      } else {
        newW <- tempW
      }
        S <- Z %*% newW
      if(PIT) {
        pS <- est.PIT(S=S,bw=bw,adjust=adjust)
        S <- pS$Fx*sqrt(12)
        fS <- pS$fx*sqrt(12)
        rm(pS)
      }
      
      newF <- multidcov(S,symmetric=symmetric)
    }
    deltaF <- curF-newF
    curF <- newF
    oldW <- newW
    rowTable <- c(iter, curF, alpha)
    if(verbose) message('iter: ',rowTable[1],'; newF: ',round(rowTable[2],6),'; alpha: ',alpha)
    Table <- rbind(Table, rowTable)
    if(symmetric) deltaS <- grad.mdCov(S,method=method)*fS else deltaS <- grad.adCov(S,method=method)*fS
    deltaW <- crossprod(Z,deltaS)
    iter <- iter+1
    if (deltaF < eps) break
  }

  convergence <- 1*(deltaF < eps)
  if(alpha < alpha.eps && convergence == 0) convergence = 2
  if (convergence==0) {
    warning("convergence not obtained in ", maxit, " iterations used.")
  } else if (convergence==2) {
    warning("check convergence: alpha is less than alpha.eps, so the norm of the gradient is greater than eps but probably sufficiently small")
  }
  if(symmetric) colnames(Table)=c('Iter','multidcov','alpha') else colnames(Table)=c('Iter','adcov','alpha') 
  
  if(whiten==TRUE) { 
    r <- list(S = Z%*%oldW, W = oldW, M = solve(oldW)%*%ginv(whiteU), f=curF, Table = Table, convergence = convergence)
  } else if(PIT) {
    Stemp = Z%*%oldW
    sdMat = diag(apply(Stemp,2,sd))
    inv.sdMat = solve(sdMat)
    oldW = oldW%*%inv.sdMat
    r <- list(S = Stemp%*%inv.sdMat, W = oldW, M = ginv(oldW), f=curF, Table = Table, convergence = convergence)
  }  else {
    r <- list(S = Z%*%oldW, W = oldW, M = ginv(oldW), f=curF, Table = Table, convergence = convergence)
  }
  r
}


