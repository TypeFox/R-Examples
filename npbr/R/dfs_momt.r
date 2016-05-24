 dfs_momt<-function(xtab, ytab, x, rho, k, ci=TRUE)
 {
  # conditions
  stopifnot(length(xtab)==length(ytab))
  if(length(rho)==1) rho<-rep(rho,length(x))
  if(length(k)==1) k<-rep(k,length(x))
  stopifnot(length(k)==length(x),length(rho)==length(x))
  
  # initialization   
    n=length(xtab)
    nx=length(x)
    res = rep(0,nx)
    lowCI_phimomt = rep(0,nx)
    upCI_phimomt = rep(0,nx)
        
   for (j in 1:nx)
    {
    # value of interest:
     xj = x[j]

    # Transformed data (z=column vector):
     z = ytab * as.numeric(xtab <= xj)
     zsort = sort(z)

    # Moments Mn1 & Mn2
     gk=zsort[n-k[j]]
     som1=0
           
     for(i in 0:(k[j]-1))
     {gi=zsort[n-i]
      som1 = som1 + (log(gi)-log(gk))
     }
     Mn1=som1/k[j]

     res[j] =  gk + (gk * Mn1 * (1 + rho[j]))
        
     var2=(rho[j]^2) * (1 + (2/rho[j]))^(-1)   
      
     corstd= gk*Mn1*(1 + 1/rho[j]) 
     lowCI_phimomt[j]= res[j] -1.96*sqrt(var2/(k[j]))*corstd
     upCI_phimomt[j] = res[j] +1.96*sqrt(var2/(k[j]))*corstd            
  }  

   if(ci)
   {res<-cbind(res, lowCI_phimomt, upCI_phimomt)}

   return(res) 
 }