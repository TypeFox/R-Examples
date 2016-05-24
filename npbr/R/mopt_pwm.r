mopt_pwm<-function(xtab, ytab, x, a=2, rho, wind.coef=0.1)
{
  # conditions
  stopifnot(length(xtab)==length(ytab))
  stopifnot(length(a)==1)
  if(length(rho)==1) rho<-rep(rho,length(x))
  stopifnot(length(rho)==length(x))
  
  # Sample and grid sizes:
  n=length(xtab)
  nx=length(x)

  # Sort the data points w.r.t the inputs Xi:
  ytab<-ytab[order(xtab)]
  xtab<-sort(xtab)

  Nx<-unlist(lapply(x,function(y) length(which(xtab<=y))))
  kxmax=unlist(lapply(floor(sqrt(Nx)),min,10))
 # kxmax=floor(Nx/4)
  nkk=kxmax[nx]
 # nkk=150
 # kxmax<-rep(nkk,nx)
  
  hatphi_tkk<-matrix(0,nkk,nx)

  # Computatuion of the desired output:
  for(j in 1:nx)
   {# value of interest:
    xj = x[j]
    evi=rho[j]

    for(kk in 1:kxmax[j])
     {
      hatphi_tkk[kk,j]=dfs_pwm(xtab, ytab, xj, kk, a, evi, ci=FALSE)
     }
   }
   
   CV=matrix(0,nkk,nx)
#  windphi_x=sapply(floor(kxmax/20),function(x) max(x,3))
  windphi_x=sapply(floor(wind.coef*kxmax),function(x) max(x,3))

  for(j in 1:nx)
  {wind=windphi_x[j]
   if((wind+1)<(kxmax[j]-wind))
    {
     for(kk in (wind+1):(kxmax[j]-wind))
      {
       CV[kk,j]=sd(hatphi_tkk[(kk-wind):(kk+wind),j])
      }
    }
  }

 CV[hatphi_tkk <= 0]= 0
 CV[CV==0]=Inf
 CVmin=min(CV)
 kopt=apply(CV,2,which.min)
 return(kopt)
}