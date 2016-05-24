ci_mean_min <- function(x, indic_col, alpha, beta)
{
  x_num   = x[,indic_col]  
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit  <- dim(as.matrix(x_num))[1]
  
  # Numeric check
  if (n_indic<2)
  {
    stop(paste("There must be at least two simple indicators!"))
  }
  if (alpha>1)
  {
    stop(paste("Alpha must be set >=0 and <=1!"))
  }
  if (alpha<0)
  {
    stop(paste("Alpha must be set >=0 and <=1!"))
  }
  if (beta<0)
  {
    stop(paste("Beta must be set >=0!"))
  }
  
  
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
  

  #punto A
  Ma <- colMeans(x_num) 
  Sqm <- (apply(x_num,2,var))^0.5
  S=10
  M=100
  z = ((x_num-Ma)/Sqm)*S + M
  
  #punto B
  Ma_z <- apply(z,1,mean)
  min_col <- apply(z,1,min)
  part = ((Ma_z - min_col)^2+beta^2)^0.5
  MMF = Ma_z - alpha * (part - beta)
  
  r<-list(ci_mean_min_est=MMF, ci_method="mean_min")
  r$call<-match.call()
  class(r)<-"CI"
  return(r)

}





