ci_wroclaw <- function(x,indic_col)
{
  options(warn=-1)
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
  
  #passo 2
  Ma <- matrix(0, nrow=1, ncol=n_indic)
  for (i in seq(1,n_indic)) 
  {
    Ma[i] <- max(x_num[,i], na.rm=TRUE)
  }
  #Ma <- apply(x_num,2,max)
  Ma_m <- t(matrix(Ma, nrow=n_indic, ncol=n_unit))
  
  #passo 3
  diff_Ma = x_num - Ma_m
  
  #passo 4
  diff_Ma_q = diff_Ma^2
  
  #passo 5
  Sum_riga <- as.matrix(apply(diff_Ma_q,1,sum))
  
  # passo 6
  Sum_riga_rad <- Sum_riga^0.5
  
  # passo 7
  media = mean(Sum_riga_rad, na.rm=TRUE)
  std   = (var(Sum_riga_rad, na.rm=TRUE))^0.5
  div <- t(matrix(media + 2*std, nrow=1, ncol=n_unit))
  ci_wroclaw_est = Sum_riga_rad / div

  r<-list(ci_wroclaw_est=ci_wroclaw_est, ci_method="wroclaw")
  r$call<-match.call()
  class(r)<-"CI"
  return(r)
  
}


