loc_max<-function(xtab,ytab,x,h,type="one-stage")
{
 stopifnot(type%in%c("one-stage","two-stage"))
 
nx = length(x)
res = matrix(0,1,nx)

 for(j in 1:nx)
 {
  ztab = ytab * as.numeric(abs(xtab-x[j])<=h)
  res[j] = max(ztab)
 }

 if(type=="two-stage")
 { for(j in 1:nx)
   {
    res[j]<-dea_est(xtab, ifelse(abs(xtab-x[j])<=h,res[j],ytab), x[j], type = "dea")
   }
 }
return(res)

}