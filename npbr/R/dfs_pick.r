 dfs_pick<-function(xtab, ytab, x, k, rho, ci=TRUE)
 {
  # conditions
  stopifnot(length(xtab)==length(ytab))
  if(length(rho)==1) rho<-rep(rho,length(x))
   stopifnot(length(rho)==length(x))
  if(length(k)==1) k<-rep(k,length(x))
  stopifnot(length(k)==length(x))
 
  # initialization
    n=length(xtab)
    nx=length(x)
    res = rep(0,nx)
    lowCI_phitilde = rep(0,nx)
    upCI_phitilde = rep(0,nx)

   for (j in 1:nx)
    {
    # value of interest:
     xj = x[j]

    # Transformed data (z=column vector):
     z = ytab * as.numeric(xtab <= xj)
     zsort = sort(z)

     kx=floor((length(which(z>0)))^0.5)- k[j] + 1

     gk1=zsort[n-kx+1]
     gk2=zsort[n-2*kx+1]

     evi=rho[j]
     var3=(evi^(-2))*(2^(-2/evi))/((2^(-1/evi) -1)^4)
     
     res[j] =  gk1 + (gk1-gk2)/(2^(1/evi)-1)
     lowCI_phitilde[j]= res[j] -1.96*sqrt(var3/(2*kx))*(gk1-gk2)
     upCI_phitilde[j] = res[j] +1.96*sqrt(var3/(2*kx))*(gk1-gk2)
  }

   if(ci)
   {res<-cbind(res,lowCI_phitilde,upCI_phitilde)}

   return(res)
 }