rho_pwm<-function(xtab, ytab, x, a=2, lrho=1, urho=Inf)
{
  # conditions
  stopifnot(length(xtab)==length(ytab))
  stopifnot(length(a)==1)

  # Sample and grid sizes:
  n=length(xtab)
  nx=length(x)
  # Initialization of the desired outputs:
  rho_x = rep(0,nx)

  # Sort the data points w.r.t the inputs Xi:
  ytab<-ytab[order(xtab)]
  xtab<-sort(xtab)
  
  pb <- txtProgressBar() 
  # Computatuion of the desired output:
  for(k in 1:nx)
   {
   # Fix the evaluation point x(k)
    xk=x[k]
    # Calculate the effective sample size N_x(k)=
    ppk=which(xtab<=xk)
    Nxk=length(ppk);
    
  # kn is the grid of possible values of coefm
  # kxmax<-floor(Nxk/4)
   kxmax<-150
 #  kxmax<-min(100,floor(Nxk/4))
   kn<-1:kxmax
   rhokn<-rep(0,length(kn))
   
   for(j in 1:kxmax)
   {
    coefm<-kn[j]
    # The order m should be defined in terms of N_x(k)
    m=max(3,round(coefm*(Nxk^(1/3))))
    # Compute the mfrontiers at x(k):
    Yxsort=sort(ytab[ppk])
    phim1<-phim2<-phim4<-max(Yxsort)
   # for(i in 1:(Nxk-1))
   #   {phim1=phim1-((i/Nxk)^m * (Yxsort[i+1]-Yxsort[i]))
   #    phim2=phim2-((i/Nxk)^(a*m) * (Yxsort[i+1]-Yxsort[i]))
   #    phim4=phim4-((i/Nxk)^((a^2)*m) * (Yxsort[i+1]-Yxsort[i]))
   #   }
    phim1<-max(Yxsort)-sum((1:(Nxk-1)/Nxk)^m*diff(Yxsort))
    phim2<-max(Yxsort)-sum((1:(Nxk-1)/Nxk)^(a*m)*diff(Yxsort))  
    phim4<-max(Yxsort)-sum((1:(Nxk-1)/Nxk)^(a*a*m)*diff(Yxsort))    
    rhokn[j]<-log(a)/log((phim1-phim2)/(phim2-phim4))
   }

    flag<-which(rhokn>lrho & rhokn<urho)
    rhokk_j<-rhokn[flag]
    kx_j=kn[flag]

    wind=max(floor(sqrt(kxmax)),3)
    nkx=length(rhokk_j)

    if(nkx==0)  # special case where no evi rhokn is superior to 1
    {rho_x[k]=NA
    }else
    {
     if(nkx < (2*wind + 1))
      {kopt_j=kx_j[nkx]
       rhohat_j=rhokk_j[nkx]
      }else
      {
      CrV=NULL
        for(kk in (wind+1):(nkx-wind))
           { CrV=c(CrV,sd(rhokk_j[(kk-wind):(kk+wind)]))}
       CVmin<-min(CrV)
       Iopt<-which.min(CrV)
       kopt_j=kx_j[Iopt+wind]
       rhohat_j=rhokk_j[Iopt+wind]
      }
    rho_x[k]=rhohat_j
   }
   setTxtProgressBar(pb, k/nx)
  }
  close(pb)
 return(rho_x)
}