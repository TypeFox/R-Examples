cub_spline_kn<-function(xtab,ytab, method,krange=1:20,type="AIC")
 {
  stopifnot(type %in% c("AIC","BIC"),length(xtab)==length(ytab), method%in%c("u","m","mc"), is.integer(krange))
       
  criteria<-NULL
 
  if(type=="AIC")
  {for (k in krange)
   {criteria<-c(criteria,log(sum(abs(ytab-cub_spline_est(xtab,ytab,xtab,k,method))))+(k+3)/length(xtab))}
  } else
  {for (k in krange)
   {criteria<-c(criteria,log(sum(abs(ytab-cub_spline_est(xtab,ytab,xtab,k,method))))+log(length(xtab))*(k+3)/(2*length(xtab)))}
  }
  
 return(krange[which.min(criteria)])
 }