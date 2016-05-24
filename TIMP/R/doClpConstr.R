"doClpConstr" <-
  function (X, clp_ind, clpCon, clpequ, num_clpequ, usecompnames0, 
            usecompnamesequ)  
  {   
    
    #add equality constraints
    if(num_clpequ > 0) {
      for(i in 1:num_clpequ) {
        X[, clpCon$clpMod[clp_ind, i] ] <- 
          X[, clpCon$clpMod[clp_ind, i] ] + 
          X[, clpCon$clpRem[clp_ind,i] ] * clpequ[i]
      }
    }
    ## column indices that are removed because of 0 constraints
    if(usecompnames0)
      rem0 <- colnames(X)[ which(clpCon$clp0mat[clp_ind, ] > 0)]
    else 
      rem0 <- which(clpCon$clp0mat[clp_ind, ] > 0)
    ## column indices that are removed because of equ constraints
    if(usecompnamesequ)
      remE <- colnames(X)[clpCon$clpRem[clp_ind, ]]
    else 
      remE <- clpCon$clpRem[clp_ind, ]
    rem <- append(rem0, remE) 
    
    if(sum(rem) > 0 )
      X <- as.matrix(X[, - rem])
    
    
    X
  }

