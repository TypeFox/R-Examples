"weightPsi" <-
  function (model) 
  {
    psisim <- model@psi.df
    if(is.na(model@weightM[1,1])) # weight not used from previous analysis 
      weight <- matrix(1, nrow=model@nt, ncol=model@nl)
    else weight <- model@weightM
    wt <- model@weightpar
    typ <- model@weightsmooth$funct
    if(length(wt)  !=  0 && length(wt$poisson) == 0) {
      for (i in 1:length(wt)) {
        if (is.na(wt[[i]][1])) 
          wt[[i]][1] <- model@x[1]
        if (is.na(wt[[i]][2])) 
          wt[[i]][2] <- model@x[model@nt]
        if (is.na(wt[[i]][3])) 
          wt[[i]][3] <- model@x2[1]
        if (is.na(wt[[i]][4])) 
          wt[[i]][4] <- model@x2[model@nl]
        for (j in 1:model@nl) {
          for (k in 1:model@nt) {
            
            if (((model@x2[j] >= wt[[i]][3]) && (model@x2[j] <= 
                                                   wt[[i]][4])) && ((model@x[k] >= wt[[i]][1]) && 
                                                                      (model@x[k] <= wt[[i]][2]))) {
              weight[k, j] <- wt[[i]][5] * weight[k, j]                    
            }
          }
        }
      }
    }
    else {
      if(wt$poisson){
        for (j in 1:model@nl) {
          for (k in 1:model@nt) { 
            if(psisim[k, j]>0){
              weight[k, j] <- 1 / (sqrt(psisim[k,j])) 
              
            }	}	 
        }
      }
    }
    psisim <- psisim * weight
    list(psi.weight = psisim, weight = weight)
  }


