 poly_degree<-function(xtab, ytab, prange=0:20, type="AIC")
 {
  stopifnot(type %in% c("AIC","BIC"),length(xtab)==length(ytab), is.integer(prange))
       
  criteria<-NULL
 
  if(type=="AIC")
  {for (deg in prange)
   {criteria<-c(criteria,log(sum(abs(ytab- poly_est(xtab,ytab,xtab,deg))))+(deg+1)/length(xtab))}
  } 
  else
  {for (deg in prange)
   {criteria<-c(criteria,log(sum(abs(ytab-poly_est(xtab,ytab,xtab,deg))))+log(length(xtab))*(deg+1)/(2*length(xtab)))}
  }
  
 return(prange[which.min(criteria)])
 }