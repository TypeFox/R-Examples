setpath <-
function(d1,d2,M=1,transform=NULL,verbose=FALSE,minalpha=NULL,normalize=TRUE,pvalue="chisq",npermutations=10000)
{
  p=dim(d1)[2]
  n1 = dim(d1)[1]
  n2 = dim(d2)[1]
  
  ## 00: normalize data to median eigenvalue if called for:
  if(normalize)
  {
    #if(p>max(n1,n2))
    #{
      e1 = eigen(cov(d1),symmetric=TRUE,only.values=TRUE)$values
      e2 = eigen(cov(d2),symmetric=TRUE,only.values=TRUE)$values
      if(n1>p){medianeigen1 = median(e1)}
      if(n2>p){medianeigen2 = median(e2)}
      if(n1<=p){medianeigen1 = median(e1[e1>1e-12])*n1/p}
      if(n2<=p){medianeigen2 = median(e2[e2>1e-12])*n2/p}
      scaling.factor = mean(medianeigen1,medianeigen2)
      d1 = d1/sqrt(scaling.factor)
      d2 = d2/sqrt(scaling.factor)
    #}
    #if(p<=max(n1,n2))
    #{	
    #  e1 = eigen(cov(d1),symmetric=TRUE,only.values=TRUE)$values
    #  e2 = eigen(cov(d2),symmetric=TRUE,only.values=TRUE)$values
    #  medianeigen1 = median(e1)
    #  medianeigen2 = median(e2)
    #  scaling.factor = mean(medianeigen1,medianeigen2)
    #  d1 = d1/sqrt(scaling.factor)
    #  d2 = d2/sqrt(scaling.factor)
    #}
  }
  ## 0: process input:
  p=dim(d1)[2]
  n1 = dim(d1)[1]
  n2 = dim(d2)[1]
  n=c(n1,n2)
  w = n/sum(n)
  d = list(d1,d2)
  
  g1 = p/n1; g2 = p/n2; g = c(g1,g2) #gammas
  covs = list()
  covs[[1]] = cov(d1)
  covs[[2]] = cov(d2)
  
  ## A: get basic stats:
  sighat = list(cov(d1),cov(d2))
  e1 = eigen(sighat[[1]],symmetric=TRUE,only.values=TRUE)$values
  e2 = eigen(sighat[[2]],symmetric=TRUE,only.values=TRUE)$values
  #L1 = e1[1]; L2 = e2[1]; L = c(L1,L2)
  L = cbind(e1[1:M],e2[1:M])
  T1 = sum(e1); T2 = sum(e2); T = c(T1,T2)
  
  ## B: estimate alpha0: the common first eigen under H0, and use it to get correction factor for L1-L2:
  # get adjusted first M eigens:
  alphabar = matrix(NA,M,2)
  a0 = QLcorrection = c()
  for(m in 1:M)
  {
    eigencorrect = unbias.eigens(L[m,],g,w,minalpha)
    alphabar[m,] = eigencorrect$a
    a0[m] = eigencorrect$a0
    QLcorrection[m] = eigencorrect$QLcorrection
  }
  
  ## B.1: estimate the number of spiked eigenvalues:
  thresh = (1 + sqrt(g))^2 + sqrt(2*log(n)/n)
  #m = max(sum(e1>thresh),sum(e2>thresh),1)
  mhat = c()
  mhat[1] = max(sum(e1>thresh[1]),M)
  mhat[2] = max(sum(e2>thresh[2]),M)   # previously: "thresh" was missing the [2] or [1].  fixed 9-7.
  
  ## B.2: get adjusted further eigens:
  spikes = list(); for(k in 1:2){spikes[[k]] = rep(NA,mhat[k])} # de-biased estimates of spiked eigenvalues
  #nullspikes = rep(NA,mhat) # estimates of common m^th eigenvalue
  #QLcorrections = rep(NA,mhat)
  
  for(k in 1:k){
    for(m in 1:mhat[k])
    {
      tempeigencorrect = unbias.eigens(c(e1[m],e2[m]),g,w,minalpha)
      spikes[[k]][m] = tempeigencorrect$a[k]
      #nullspikes[i] = tempeigencorrect$a0
    }}
  
  ## C: estimate Sigma_Q:
  covQ = matrix(0,M+1,M+1)
  
  ## Ci: get var(T1), var(T2) using theory and assuming H0 
  # (actually, assuming all vars are independent for convenience now.  fix later)
  varT = c()
  for(k in 1:2)
  {
    varT[k] =  2*(sum(a0^2)/n[k] + (p-M)/n[k])
    if(mhat[k]>M){ varT[k] = varT[k] + 2/n[k]*(sum( (spikes[[k]][(M+1):mhat[k]])^2) - (mhat[k]-M) )}  # if additional spikes, adjust the above result appropriately.
  }
  covQ[M+1,M+1]=sum(varT)
  
  ## Cii: get var(L1-L2-bL):
  for(m in 1:M)
  {
    rho = theta = varLs = c()
    c0 = (1/g[1]+1/g[2])^2*(a0[m]-1)^2
    for(k in 1:2)
    {
      rho[k] = a0[m]*(1 + g[k]/(a0[m]-1))
      deriv.f.k = .5*(1 + (rho[k] - 1 - g[k])/sqrt((rho[k]-1-g[k])^2-4*g[k]))
      theta[k] = 1 + (g[1]-g[2])/c0*deriv.f.k/g[k]
      varLs[k] = 2*a0[m]/n[k]*theta[k]^2*rho[k]/(1+a0[m]*g[k]/((a0[m]-1)^2-g[k]))
    }
    covQ[m,m] = sum(varLs)
  }
  
  ## Ciii: get cov(L,T)   
  for(m in 1:M)
  {
    rho = theta = covLTs = c()
    c0 = (1/g[1]+1/g[2])^2*(a0[m]-1)^2
    for(k in 1:2)
    {
      rho[k] = a0[m]*(1 + g[k]/(a0[m]-1))
      deriv.f.k = .5*(1 + (rho[k] - 1 - g[k])/sqrt((rho[k]-1-g[k])^2-4*g[k]))
      theta[k] = 1 + (g[1]-g[2])/c0*deriv.f.k/g[k]
      covLTs[k] = 2*a0[m]/n[k]*theta[k]*rho[k]/(1+a0[m]*g[k]/((a0[m]-1)^2-g[k]))
    }
    covLT = sum(covLTs) 
    covQ[m,M+1]=covQ[M+1,m]=covLT
  }
  
  ## D: compile Q:
  Q = c(L[,1] - L[,2] - QLcorrection, T1-T2)
  
  ## D: return test stat (chisq^2_2 under H0):
  # first: linear transformation: call it I if no argument was entered.
  A = transform
  if(length(transform)==0){A = diag(M+1)}
  A=as.matrix(A)
  if(dim(A)[1]!=M+1){stop("Dimension of linear transformation does not match the dimension of the null hypothesis.")}
  stat = t(Q) %*% A %*% solve(t(A) %*% covQ %*% A) %*% t(A) %*% Q
  out = list()
  out$stat = stat
  if(pvalue=="chisq"){	out$pval = 1-pchisq(stat,dim(A)[2]) }
  
  if(pvalue=="permutation")
  {
    d.combined = rbind(d1,d2)
    permstats = c()
    for(i in 1:npermutations)
    {
      prows1 = sample(1:dim(d.combined)[1],dim(d1)[1],replace=FALSE)
      prows2 = setdiff(1:dim(d.combined)[1],prows1)
      permstats[i] = setpath(d1=d.combined[prows1,],d2=d.combined[prows2,],M=M,transform=transform,verbose=FALSE,minalpha=minalpha,normalize=FALSE,pvalue="chisq")$stat			
    }
    out$pval = mean(as.vector(stat)<permstats)
  }
  
  if(verbose)
  {
    out$stats = rbind(L,c(T1,T2))
    out$a0 = a0
    out$correction = QLcorrection
    out$covQ = covQ
    out$m = m
  }
  return(out)
}
