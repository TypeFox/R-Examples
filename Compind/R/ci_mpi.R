ci_mpi <- function(x, indic_col, penalty="POS")
{
  x_num   = x[,indic_col]  
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit  <- dim(as.matrix(x_num))[1]
  
  # Numeric check
  if (n_indic<2)
  {
    stop(paste("There must be at least two simple indicators!"))
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
  Sqm_z <- (apply(z,1,var))^0.5
  cv = Sqm_z / Ma_z
  
  #punto C
  if (penalty=="POS") {
    ci_mpi_est <- Ma_z*(1-cv^2)
  } else {
    ci_mpi_est <- Ma_z*(1+cv^2)
  }

  r<-list(ci_mpi_est=ci_mpi_est, ci_method="mpi")
  r$call<-match.call()
  class(r)<-"CI"
  return(r)

  
}





