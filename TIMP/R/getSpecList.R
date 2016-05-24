"getSpecList" <- 
  function (multimodel, t, getclperr = FALSE) 
  {   
    m <- multimodel@modellist
    clpList <- list()
    resultlist <- multimodel@fit@resultlist
    for (i in 1:length(m)) {
      ## fill in the clp 
      if(getclperr) {
        cptemp <- resultlist[[i]]@std_err_clp
      } else {
        cptemp <- resultlist[[i]]@cp
      }
      
      if(m[[i]]@clpType == "x2") {
        nx <- m[[i]]@nl
        colc <- max(m[[i]]@ncolc)
      }	
      else {
        nx <- m[[i]]@nt
        colc <- max(m[[i]]@ncole)
      }
      X <- matrix(nrow = nx, ncol = colc)
      
      X[which(m[[i]]@clpCon$clp0mat!=0)] <- 0
      
      ## mark those clp determined by equality with 0 
      if(m[[i]]@lclpequ) {
        for( k in 1:dim(m[[i]]@clpCon$clpRem)[1]) {
          for( l in 1:dim(m[[i]]@clpCon$clpRem)[2]) {
            if(m[[i]]@clpCon$clpRem[k,l]!=0) {
              X[k,m[[i]]@clpCon$clpRem[k,l]] <- 0
            }
          }
        }
      }
      
      
      for(j in 1:nx) {
        X[j,is.na(X[j,])] <- cptemp[[j]]
        if(m[[i]]@lclpequ) {
          for( l in 1:dim(m[[i]]@clpCon$clpRem)[2]) {
            if(m[[i]]@clpCon$clpRem[j,l]!=0) {
              X[j,m[[i]]@clpCon$clpRem[j,l]] <- X[j,m[[i]]@clpCon$clpMod[j,l]] * t[[i]]@clpequ[l]
            }
          }
        }
      }
      
      ## go back and enforce the equ constraints
      # this is now done above:
      #if(m[[i]]@lclpequ) {
      #for(j in 1:nx) {
      #cptemp[[k]][m[[i]]@clpCon$clpMod[k,l]] * t[[i]]@clpequ
      #  X[j, m[[i]]@clpCon$clpRem[j,] ] <- X[j, m[[i]]@clpCon$clpMod[j,] ] *
      #    t[[i]]@clpequ
      #}
      clpList[[i]] <- X
    }  
    clpList
  }
