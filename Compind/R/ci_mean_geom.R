ci_mean_geom <- function(x, indic_col, na.rm=TRUE)
{
 # library(psych)
  x_num   = x[,indic_col]
  n_indic <- dim(x_num)[2]
  n_unit  <- dim(x_num)[1]
  
  # Numeric check
  for (i in seq(1,n_indic)) 
  {
    if (!is.numeric(x_num[,i]))
    {
      stop(paste("Data set not numeric at column:",i))
    }
  }   

  
  for (i in seq(1,n_unit)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (is.na(x_num[i,j]))
      {
        message(paste("Pay attention: NA values at column:",i,", row",j,". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
        #       options(warn=-2)  
      }
    }
  }   
  
    
  
  t_x = t(x_num)
  
  if (!na.rm) {
  ci_mean_geom_est <- geometric.mean(t_x)
  } else {
  ci_mean_geom_est <- geometric.mean(t_x,na.rm=FALSE)  
  }
  
  r<-list(ci_mean_geom_est=ci_mean_geom_est, ci_method="mean_geom")
  r$call<-match.call()
  class(r)<-"CI"
  return(r)
  
}

