ci_bod_dir <- function(x,indic_col,dir)
{
#   library(Benchmarking)
#   library(MASS)
#   library(Hmisc)
  
  x_num   = x[,indic_col]
  ci_data = as.matrix(cbind(x_num))
  ci_data_dir= as.matrix(cbind(x_num))
  n_indic = dim(x_num)[2]
  n_unit  = dim(x_num)[1]
  uni <- as.matrix(seq(1, 1, len = n_unit))
  
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
  

  
  for (i in seq(1,n_indic)) 
       {
           ci_data_dir[,i] = ci_data[,i] * dir[i]
       }
    
  ci_dir  = Benchmarking::dea(uni,ci_data_dir,RTS="crs", ORIENTATION="out", DIRECT=dir)
  ci_bod_dir_est = 1/(ci_dir$eff + 1)
  
  r<-list(ci_bod_dir_est=ci_bod_dir_est, ci_method="bod_dir")
  r$call<-match.call()
  class(r)<-"CI"
  r
  
}

