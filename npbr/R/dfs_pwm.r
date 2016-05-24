dfs_pwm<-function(xtab, ytab, x, coefm, a=2, rho, ci=TRUE)
{
  # conditions
  stopifnot(length(xtab)==length(ytab))
  if(length(rho)==1) rho<-rep(rho,length(x))
  if(length(coefm)==1) coefm<-rep(coefm,length(x))
  stopifnot(length(a)==1,length(rho)==length(x), length(coefm)==length(x)) 

  # Sample and grid sizes: 
  n=length(xtab)
  nx=length(x)
  # Initialization of the desired outputs: 
  res = rep(0,nx)
  Nxk <- rep(0, nx)
  m <- rep(0, nx)
  lowCI_res = rep(0,nx)
  upCI_res = rep(0,nx)
  YXSORT<-NULL       
  # Sort the data points w.r.t the inputs Xi:
  ytab<-ytab[order(xtab)]
  xtab<-sort(xtab)

  # Computatuion of the desired outputs:
  for(k in 1:nx)
   {
    # Fix the evaluation point x(k):
    xk=x[k]
    # Calculate the effective sample size N_x(k)=
    ppk=which(xtab<=xk)
    Nxk[k]=length(ppk)
    # The order m should be defined in terms of N_x(k)
    m[k]=max(3,round(coefm[k]*(Nxk[k]^(1/3))))
    # Compute the mfrontiers at x(k):
    Yxsort=sort(ytab[ppk])

     phim1<-max(Yxsort)-sum((1:(Nxk[k]-1)/Nxk[k])^m[k]*diff(Yxsort))
     phim2<-max(Yxsort)-sum((1:(Nxk[k]-1)/Nxk[k])^(a*m[k])*diff(Yxsort))

    # Fix the extreme-value index rho(k):
    rhok=rho[k]
    # compute hat{ell}_xk and tilde{phi}_m(xk):
    aatr=m[k]^(-1/rhok)
    abtr=gamma(1 + 1/rhok)
    actr=(1-(1/a)^(1/rhok))
    # hat{ell}_xk:
    ell=( (aatr*abtr*actr)/(phim2-phim1) )^(rhok)
    # tilde{phi}_m(xk):
    Bias = abtr*(1/(m[k]*ell))^(1/rhok)
    res[k] = phim1 + Bias 
    
    YXSORT<-c(YXSORT,Yxsort)
   }

  if(ci)
  {# Empirical asymptotic variance \hat\sigma^2(m,x):
    sigma2= .C("sigma2m",as.double(Nxk), as.double(m), as.integer(n), as.integer(nx), as.double(YXSORT), sigma2=as.double(rep(0,nx)),PACKAGE="npbr") 
    sigma2=sigma2$sigma2

   # 95% confidence interval of the frontier:
    lowCI_res= res - 1.96*sqrt(sigma2/n)
    upCI_res = res + 1.96*sqrt(sigma2/n)   
    res<-cbind(res, lowCI_res, upCI_res)
   }

   return(res) 
}