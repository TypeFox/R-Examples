ci_rbod <- function(x,indic_col,M,B)
{
 # library(nonparaeff)
  x_num   = x[,indic_col]
  n_indic = dim(x_num)[2]
  n_unit  = dim(x_num)[1]
  uni <- as.matrix(seq(1, 1, len = n_unit))
  data_per_Rbod = cbind(x_num,uni)
  
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
  
  
  
  stima_m <- nonparaeff::orderm(data_per_Rbod, noutput = n_indic,
                                orientation = 2, M=M, B=B)
  ci_rbod_est = 1/stima_m[[1]]
  
  r<-list(ci_rbod_est=ci_rbod_est, ci_method="rbod")
  r$call<-match.call()
  class(r)<-"CI"
  return(r)
  
}
