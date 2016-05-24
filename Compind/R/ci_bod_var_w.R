ci_bod_var_w <- function(x,indic_col,boot_rep = 5000)
{
#   library(Hmisc)
#   library(lpSolve)
#   library(nonparaeff)
#   library(boot)
  
  boot_r = boot_rep
  x_num   = x[,indic_col]
  n_indic = dim(x_num)[2]
  n_unit  = dim(x_num)[1]
  
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
  
  
  # lower and higher variance confidence limits #################
  low_var  = matrix(0, nrow=n_indic, ncol=1)
  high_var = matrix(0, nrow=n_indic, ncol=1)
  for (i in seq(1,n_indic)) 
  {
    bootcorr <- boot(x_num[,i], var, R=boot_r)
    interv_var = boot.ci(bootcorr, conf = 0.95, type="norm")
    low_var[i]  = interv_var$normal[2]
    high_var[i] = interv_var$normal[3]
  }
  
  for (i in seq(1,n_indic)) 
  {
    if (low_var[i]<=0) {
      stop("Warning: Lower confidence bound negative!") 
    } 
  }
  
  # lower and higher ratio matrices ############################
  conf_low = matrix(0, nrow=n_indic, ncol=n_indic)
  for (i in seq(1,n_indic)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (i<j)
      { 
        conf_low[i,j] = low_var[i] / high_var[j]
      }
    }
  }
  
  conf_high = matrix(0, nrow=n_indic, ncol=n_indic)
  for (i in seq(1,n_indic)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (i<j)
      { 
        conf_high[i,j] = high_var[i] / low_var[j]
      }
    }
  }
  
  ## Caricamento matrice restriction weights ################
  
  num_righe_weights = ((n_indic*n_indic)-n_indic)/2
  ## (LEFT lower bounds)
  left_low = matrix(0, nrow=num_righe_weights, ncol=n_indic)
  riga = 0
  for (i in seq(1,n_indic)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (conf_low[i,j]>0)
      { 
        riga = riga + 1
        left_low[riga,i] =1
        left_low[riga,j] =-conf_low[i,j]
      }
    }
  }
  
  ## (LEFT higher bounds)
  left_high = matrix(0, nrow=num_righe_weights, ncol=n_indic)
  riga = 0
  for (i in seq(1,n_indic)) 
  {
    for (j in seq(1,n_indic)) 
    {
      if (conf_high[i,j]>0)
      { 
        riga = riga + 1
        left_high[riga,i] =-1
        left_high[riga,j] =conf_high[i,j]
      }
    }
  }
  
  # LEFT - Colonna di zero per l'input unitario
  left_zero = matrix(0, nrow=num_righe_weights*2, ncol=1)
  
  #LEFT - Matrice complessiva
  left_low_high = rbind(left_low,left_high)
  left          = cbind(left_low_high,left_zero)
  
  #RIGHT - condizioni di zero
  right = matrix(0, nrow=num_righe_weights*2, ncol=1)
  
  #CONDIZIONE
  dir = matrix(">=", nrow=num_righe_weights*2, ncol=1)
  
  ########  creazione input unitario
  uni <- as.matrix(seq(1, 1, len = n_unit))
  data_per_bod = cbind(x_num,uni)
  
  CI = ar.dual.dea(data_per_bod, noutput = n_indic, orientation = 2, rts = 1, 
                   ar.l = left, 
                   ar.r = right,
                   ar.dir = dir)
  
  ci_bod_var_w_est = 1/CI[[1]]
  
  r<-list(ci_bod_var_w_est=ci_bod_var_w_est, ci_method="bod_var_w")
  r$call<-match.call()
  class(r)<-"CI"
  r
  
}

