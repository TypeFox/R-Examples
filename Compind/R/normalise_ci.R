normalise_ci <- function(x, indic_col, polarity, method=1, z.mean=0, z.std=1, ties.method ="average")
  {
  
  ci_norm_finale = x
  x_num   = x[,indic_col]
  n_indic <- dim(x_num)[2]

  # Numeric check
  for (i in seq(1,n_indic)) 
    {
      if (!is.numeric(x_num[,i]))
        {
          stop(paste("Data set not numeric at column:",i))
        }
    }  
  
  
  # 1 condizione di errore: dati negativi
  for (i in seq(1,n_indic)) 
  {
    if (min(x_num[,i],na.rm = TRUE)<0) 
    {
      stop("Error: simple indicator must be positive!") 
    }
  }
                               
  #1 - Standardisation or z-scores ############################
  
  if (method==1) {
      ci_norm = x_num  
      Ma  <- colMeans(x_num, na.rm = TRUE)
      Sqm <- matrix(0, nrow=1, ncol=n_indic)
      for (i in seq(1,n_indic)) 
      {
        Sqm[i] <- var(x_num[,i], na.rm=TRUE)
      }
        
      for (i in seq(1,n_indic)) 
      {
        if (polarity[i]=="POS") 
        {
          ci_norm[,i] = z.mean + ((x_num[,i]-Ma[i])/Sqm[i])*z.std
        }
        if (polarity[i]=="NEG") 
        {
          ci_norm[,i] = z.mean - ((x_num[,i]-Ma[i])/Sqm[i])*z.std
        }
        if (polarity[i]!="NEG" & polarity[i]!="POS")
        {
          stop("Please check polarity!") 
        }   
      }    
  }  
  
  #2 - min-max method #########################################
  
  if (method==2) {
    ci_norm = x_num
    
     min <- matrix(0, nrow=1, ncol=n_indic)    
     for (i in seq(1,n_indic)) 
     {
       min[i] <- min(x_num[,i],na.rm=TRUE)
     }
     max <- matrix(0, nrow=1, ncol=n_indic)    
     for (i in seq(1,n_indic)) 
     {
       max[i] <- max(x_num[,i],na.rm=TRUE)
     }    
  
    for (i in seq(1,n_indic)) 
      {
        if (polarity[i]=="POS") 
          {
            ci_norm[,i] = (x_num[,i]-min[i])/(max[i] - min[i])  
        }
        if (polarity[i]=="NEG") 
        {
            ci_norm[,i] = ((max[i]) - (x_num[,i])) / ((max[i]) - (min[i]))  
        }     
        if (polarity[i]!="NEG" & polarity[i]!="POS")
        {
            stop("Please check polarity!") 
        } 
      }
  }  
  
  #3 - ranking ################################################
  
  if (method==3) {
    ci_norm = x_num
    for (i in seq(1,n_indic)) 
    {
      if (polarity[i]=="POS") 
      {
        ci_norm[,i] <- rank(x_num[,i], na.last="keep")
      }
      if (polarity[i]=="NEG") 
      {
        ci_norm[,i] <- rank(-x_num[,i], na.last="keep")
      }
      if (polarity[i]!="NEG" & polarity[i]!="POS")
      {
        stop("Please check polarity!") 
      }   
    }
  }  
  
  r<-list(ci_norm=ci_norm,norm_method=method)
  #r$call<-match.call()
  #class(r)<-"norm_CI"
  return(r)
}
