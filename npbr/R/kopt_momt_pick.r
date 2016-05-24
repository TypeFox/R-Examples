kopt_momt_pick<-function(xtab, ytab, x, rho, method="moment", wind.coef=0.1)
{
  # conditions on data
  stopifnot(method%in%c("moment","pickands"))
  stopifnot(length(xtab)==length(ytab))
  if(length(rho)==1) rho<-rep(rho,length(x))
  stopifnot(length(rho)==length(x))
   
  # initialization 
  n=length(xtab)
  nx=length(x)

  Nx<-unlist(lapply(x,function(y) length(which(xtab<=y))))
  kxmax=floor((Nx)^(1/2))
  nkk=kxmax[nx]
  
  hatphi_tkk<-matrix(0,nkk,nx)
 
  if(method=="pickands")
   {for(j in 1:nx)     
    {# value of interest:
     xj = x[j]
     # Transformed data (z=column vector):
     z = ytab * as.numeric(xtab <= xj)
     zsort = sort(z)
   
     for(kk in 1:kxmax[j]-1)
     { 
      kx=kxmax[j] - kk + 1                        
      if(kx < 2) {break}        
      gk1=zsort[n-kx+1]
      gk2=zsort[n-2*kx+1]
      denom=gk1-gk2      
      evi=rho[j]               
      hatphi_tkk[kk,j]=gk1 + (gk1-gk2)/(2^(1/evi)-1)                  
     }  
    }   
   } 
  else
  {for(j in 1:nx)     
   {# value of interest:
    xj = x[j]
    evi=rho[j]  
    # Transformed data (z=column vector):
    z = ytab * as.numeric(xtab <= xj)
    zsort = sort(z)
   
    for(kk in 1:kxmax[j]-1)
     {                   
      hatphi_tkk[kk,j]=dfs_momt(xtab, ytab, xj, evi, kk, ci=FALSE)
     }                               
   }  
  }
             
  CV=matrix(0,nkk,nx)
#  windphi_x=sapply(floor(kxmax/20),function(x) max(x,3))
  windphi_x=sapply(floor(wind.coef*kxmax/2),function(x) max(x,3))

 
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