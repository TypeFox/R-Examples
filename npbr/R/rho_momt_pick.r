rho_momt_pick<-function(xtab, ytab, x, method="moment", lrho=1, urho=Inf)
{
 # conditions on data
 stopifnot(length(xtab)==length(ytab))
 stopifnot(method%in%c("moment","pickands"))
 
 # initialization
 n=length(xtab)
 nx=length(x)
 rho_x<-rep(0,nx)
 Nx<-unlist(lapply(x,function(y) length(which(xtab<=y))))

 for(j in 1:nx)
  {
   # value of interest:
    xj = x[j]
   # Transformed data (z=column vector):
    z = ytab * as.numeric(xtab <= xj)
    zsort = sort(z)

   # kn is the grid of possible values of k
   kxmax<-ifelse(method=="moment", ceiling(Nx[j]-1), floor(Nx[j]/4))
   kn<-1:kxmax

   # rhokn are the candidates for rho depending on kn
    rhokn<-rep(0,length(kn))

   # in this loop, we compute the values of rhokn
    if(method=="moment")
     {for(k in 1:length(kn))
      {
       gk=zsort[n-kn[k]]
       Mn1=mean(log(zsort[n:(n-kn[k]+1)])-log(gk))
       Mn2=mean((log(zsort[n:(n-kn[k]+1)])-log(gk))^2)

      # Moment's tail-index estimator:
       c1 = 1 - (Mn1^2/Mn2)
       gamman = Mn1 + 1 - (0.5/c1)
       rhokn[k]<- -(1/gamman)
      }
    }
    else
     {for(k in 1:length(kn))
      {
        gk1=zsort[n-kn[k]+1]
        gk2=zsort[n-2*kn[k]+1]
        gk4=zsort[n-4*kn[k]+1]

        num= gk2 - gk4
        denom= gk1 - gk2

        rhokn[k]<- log(2)/(log(num/denom))
       }
     }
    
    flag<-which(rhokn>lrho & rhokn<urho)
    rhokk_j<-rhokn[flag]
    kx_j=kn[flag]

    wind=max(floor(sqrt(kxmax)),3)
    nkx=length(rhokk_j)

    if(nkx==0)  # special case where no evi rhokn is superior to 1
    {rho_x[j]=NA
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
    rho_x[j]=rhohat_j
   }
  }
  return(rho_x)
}