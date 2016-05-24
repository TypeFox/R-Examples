"specparF" <- 
  function (specpar, xi, i, specref, specdispindex, specdisppar, 
            parmufunc = "") 
  {
    disppar <- specpar
    for (j in 1:length(specdispindex)) {
      if (parmufunc == "exp") {
        if (length(specdisppar[[j]]) == 1) 
          specdisppar[[j]][2] <- specdisppar[[1]][2]
        disppar[[specdispindex[[j]][1]]][specdispindex[[j]][2]] <- simpExp(specpar[[specdispindex[[j]][1]]][specdispindex[[j]][2]], 
                                                                           specdisppar[[j]], xi, specref)
      }
      if (parmufunc == "multiexp") {
        disppar[[specdispindex[[j]][1]]][specdispindex[[j]][2]] <- 
          simpExp(specpar[[specdispindex[[j]][1]]][specdispindex[[j]][2]], 
                  specdisppar[[j]], xi, specref)
        
        
      }
      if (parmufunc == "" || parmufunc == "poly") {
        disppar[[specdispindex[[j]][1]]][specdispindex[[j]][2]] <- simpPol(specpar[[specdispindex[[j]][1]]][specdispindex[[j]][2]], 
                                                                           specdisppar[[j]], xi, specref)
      }
    }
    disppar
  }
