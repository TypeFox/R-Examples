pick_est<-function(xtab, ytab, x, h, k, type="one-stage")
{
 stopifnot(length(xtab)==length(ytab),length(x)==length(k),type%in%c("one-stage","two-stage"))

n=length(xtab)
nx = length(x)
res = matrix(0,1,nx)

 for(j in 1:nx)
 {
  ztab = ytab * as.numeric(abs(xtab-x[j])<=h)
  zsort = sort(ztab)
  zk=zsort[n-k[j]]
  z2k=zsort[n-2*k[j]]
  z4k=zsort[n-4*k[j]]
  res[j] = zk + (zk-z2k)*(2^(-log((zk-z2k)/(z2k-z4k))/log(2))-1)^(-1)
 }

 if(type=="two-stage")
 { for(j in 1:nx)
   {
    res[j]<-dea_est(xtab, ifelse(abs(xtab-x[j])<=h,res[j],ytab), x[j], type = "dea")
   }
 }
 
return(res)

}