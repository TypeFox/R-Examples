########################################################################
#Function Effects
########################################################################


Effects<-function(response, factors,  effect = c("unweighted", 
                                                 "weighted")) 
{
  
  effect <- match.arg(effect)
  factorx <- factors
  samples <- split(response, factorx)
  fl <- levels(factorx)
  a <- nlevels(factorx)
  n <- sapply(samples,length)
  N <- sum(n)
  tmp1 <- sort(rep(1:a, a))
  tmp2 <- rep(1:a, a)
  
  p <- sapply(1:(a^2), function(arg) {
    x1 <- samples[[tmp1[arg]]]
    x2 <- samples[[tmp2[arg]]]
    rx1x2 <- rank(c(x1, x2))
    l1 <- length(x1)
    l2 <- length(x2)
    1/(l1 + l2) * (mean(rx1x2[(l1 + 1):(l1 + l2)]) - mean(rx1x2[1:l1])) + 
      0.5
  })
  
  intRanks <- lapply(samples, rank)
  
  V <- rep(0, a^4)
  help <- expand.grid(1:a, 1:a, 1:a, 1:a)
  h1 <- help[, 4]
  h2 <- help[, 3]
  h3 <- help[, 2]
  h4 <- help[, 1]
  for (u in 1:(a^4)) {
    i <- h1[u]
    j <- h2[u]
    r <- h3[u]
    s <- h4[u]
    if (i == r && j == s && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      ni <- length(xi)
      nj <- length(xj)
      ri <- rank(xi)
      rj <- rank(xj)
      rij <- rank(c(xi, xj))
      pj <- 1/ni * (rij[(ni + 1):(ni + nj)] - rj)
      pi <- 1/nj * (rij[1:ni] - ri)
      vi <- var(pi)/ni
      vj <- var(pj)/nj
      V[u] <- N * (vi + vj)
    }
    if (i == s && j == r && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      ni <- length(xi)
      nj <- length(xj)
      ri <- rank(xi)
      rj <- rank(xj)
      rij <- rank(c(xi, xj))
      pj <- 1/ni * (rij[(ni + 1):(ni + nj)] - rj)
      pi <- 1/nj * (rij[1:ni] - ri)
      vi <- var(pi)/ni
      vj <- var(pj)/nj
      V[u] <- -N * (vi + vj)
    }
    if (i == r && j != s && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      xs <- samples[[s]]
      ni <- length(xi)
      nj <- length(xj)
      ns <- length(xs)
      ri <- rank(xi)
      rj <- rank(xj)
      rs <- rank(xs)
      rij <- rank(c(xi, xj))
      ris <- rank(c(xi, xs))
      pij <- 1/nj * (rij[1:ni] - ri)
      pis <- 1/ns * (ris[1:ni] - ri)
      V[u] <- N * (cov(pij, pis)/ni)
    }
    if (i != r && j == s && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      xr <- samples[[r]]
      ni <- length(xi)
      nj <- length(xj)
      nr <- length(xr)
      ri <- rank(xi)
      rj <- rank(xj)
      rr <- rank(xr)
      rji <- rank(c(xj, xi))
      rjr <- rank(c(xj, xr))
      pji <- 1/ni * (rji[1:nj] - rj)
      prj <- 1/nr * (rjr[1:nj] - rj)
      V[u] <- N * (cov(pji, prj)/nj)
    }
    if (i == s && j != r && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      xr <- samples[[r]]
      ni <- length(xi)
      nj <- length(xj)
      nr <- length(xr)
      ri <- rank(xi)
      rj <- rank(xj)
      rr <- rank(xr)
      rij <- rank(c(xi, xj))
      rir <- rank(c(xi, xr))
      pij <- 1/nj * (rij[1:ni] - ri)
      pir <- 1/nr * (rir[1:ni] - ri)
      V[u] <- -N * (cov(pij, pir)/ni)
    }
    if (i != s && j == r && i != j && r != s) {
      xi <- samples[[i]]
      xj <- samples[[j]]
      xs <- samples[[s]]
      ni <- length(xi)
      nj <- length(xj)
      ns <- length(xs)
      ri <- rank(xi)
      rj <- rank(xj)
      rs <- rank(xs)
      rji <- rank(c(xj, xi))
      rjs <- rank(c(xj, xs))
      pji <- 1/ni * (rji[1:nj] - rj)
      pjs <- 1/ns * (rjs[1:nj] - rj)
      V[u] <- -N * (cov(pji, pjs)/nj)
    }
  }
  V1 <- matrix(V, ncol = a^2, nrow = a^2)
  switch(effect, weighted = {
    W <- kronecker(t(n/N), diag(a))
    samplesR <- split(rank(response), factorx)
    varsF <- sapply(samplesR,var)
    VH0F = diag(varsF/(n*N^2))
    Si2 = unlist(lapply(1:a,function(arg) var(samplesR[[arg]]-intRanks[[arg]])))
    dfBF = (sum(Si2/(N-n)))^2 / sum((Si2/(N-n))^2/(n-1)) #(sum(Si2[[arg]]/N-n[arg]))^2 / (sum((Si2[[arg]]/(N-n[arg]))^2 / (n[arg]-1))))
    text.output.W <- paste("Global Ranks")
  }, 
  unweighted = {
    W <- kronecker(t(rep(1/a, a)), diag(a))
    samplesR = lapply(1:a,function(arg){
      helpmat=rbind(1:a,matrix(1:a,nrow=a,ncol=a))
      x1 <- samples[[helpmat[1,arg]]]
      help=0
      for(j in 1:a){
        x2 <- samples[[helpmat[j+1,arg]]]
        help=help+1/length(x2)*(rank(c(x1,x2))[1:length(x1)] - rank(x1))
      }
      N/a*help+1/2})
    
    varsF <- sapply(samplesR,var)
    VH0F = diag(varsF/(n*N^2))
    Si2 = unlist(lapply(1:a,function(arg) var(samplesR[[arg]]-intRanks[[arg]])))
    dfBF = (sum(Si2/(N-n)))^2 / sum((Si2/(N-n))^2/(n-1)) #(sum(Si2[[arg]]/N-n[arg]))^2 / (sum((Si2[[arg]]/(N-n[arg]))^2 / (n[arg]-1))))
    text.output.W <- paste("Global Pseudo Ranks")
    
    
  })
  
  pd = W%*%p
  VV=W%*%V1%*%t(W)
  
  
  result <- list(pd=pd, VH0F=VH0F,VBF=VV, N=N, n=n,dfATS=dfBF)
  return(result)
}

###################################################################
# Function for Confidence Limits
###################################################################

Limits <- function(p,V,alpha,N){
  
  
  L <- p - qnorm(1-alpha/2)/sqrt(N)*sqrt(c(diag(V)))
  U <- p + qnorm(1-alpha/2)/sqrt(N)*sqrt(c(diag(V)))
  
  Psi <- diag(1/(p*(1-p)))
  VLogit <- Psi%*%V%*%t(Psi)
  Llogit <- expit(logit(p)- qnorm(1-alpha/2)/sqrt(N)*sqrt(c(diag(VLogit))))
  Ulogit <- expit(logit(p)+ qnorm(1-alpha/2)/sqrt(N)*sqrt(c(diag(VLogit))))
  
  res=list(Normal=cbind(L,U), Logit= cbind(Llogit,Ulogit))
  return(res)
}

#################################################################
#Wald Type Statistics
#################################################################
Wald <- function(M,H,V){
  WTS = t(H%*%M)%*%ginv(H%*%V%*%t(H))%*%H%*%M
  dfWTS = rankMatrix(H)[1]
  pv.WTS = 1-pchisq(WTS,dfWTS)
  res.WTS = c(WTS,dfWTS,pv.WTS)
  return(res.WTS)
}

#################################################################
#ANOVA Type Statistics
#################################################################

ANOVATYP <- function(M,H,V,n){
  
  C <- t(H)%*%ginv(H%*%t(H))%*%H
  spur <- sum(diag(C%*%V))
  D <- diag(C)*diag(ncol(C))
  Lambda <- diag(1/(n-1))
  ATS <- 1/spur*t(M)%*%C%*%M
  df_ATS1 <- spur^2/sum(diag(C%*%V%*%C%*%V))
  df_ATS2 <- spur^2/sum(diag(D%*%D%*%V%*%V%*%Lambda))
  
  pv.ATS <- 1-pf(ATS, df_ATS1, df_ATS2)
  res <- c(ATS, df_ATS1, df_ATS2, pv.ATS)
  return(res)
}

ANOVATYPH0P <- function(M,H,V,n,df){
  C <- t(H)%*%ginv(H%*%t(H))%*%H
  spur <- sum(diag(C%*%V))
  ATS <- sum(n)/spur*t(M)%*%C%*%M
  df_ATS1 <- spur^2/sum(diag(C%*%V%*%C%*%V))
  pv.ATS <- 1-pf(ATS, df_ATS1, df)
  res <- c(ATS, df_ATS1, df, pv.ATS)
  return(res)
}

#################################################################
#Logit Transformation
#################################################################

logit <- function(p){
  return(log(p/(1-p)))}
expit<-function(p){return(exp(p)/(1+exp(p)))}

################################################################
# for twosample problems
################################################################

wilcoxontests <- function(x,y){
  n1 = length(x)
  n2 = length(y)
  N=n1+n2
  
  rxy = rank(c(x,y))
  rx=rxy[1:n1]
  ry=rxy[(n1+1):N]
  W=sqrt(n1*n2/N)*(mean(ry)-mean(rx))/sd(rxy)
  
  erg=data.frame(W=W, pV2=min(2*pnorm(W),2-2*pnorm(W)),pV.L=pnorm(W),pV.U=1-pnorm(W))
  erg
}


BMstat = function(x,y,nx,ny){
  Nxy = nx + ny
  xy = c(x,y)
  rxy = rank(xy)
  rx = rank(x)
  ry= rank(y)
  plx = 1/ny*(rxy[1:nx]-rx)
  ply = 1/nx*(rxy[(nx+1):(Nxy)] -ry)
  vx = Nxy*(var(plx)/nx + var(ply)/ny)
  pd=mean(ply)
  pd0=(pd==0)
  pd1 = (pd==1)
  pd[pd0] = (1/(2*nx))/ny
  pd[pd1] = (nx-1/(2*ny))/nx
  vx0= (vx==0)
  vx[vx0] = Nxy/(2*nx^2*ny^2)
  T = sqrt(Nxy)*(pd-1/2)/sqrt(vx)
  
  #-------------Logit----------------#
  slogit=(1/(pd*(1-pd))*sqrt(vx))
  Tlogit = sqrt(Nxy)*logit(pd)/(slogit)
  
  #------------Probit-----------------#
  sprobit=sqrt(2*pi)/(exp(-1/2*(qnorm(pd))^2))*sqrt(vx)
  Tprobit = sqrt(Nxy) * qnorm(pd)/(sprobit)
  
  erg=data.frame(T=T, Logit=Tlogit,Probit=Tprobit, sdx=vx, slogit=slogit, sprobit=sprobit)
  erg
}



