z.nom.diag <-
function(Z){
    # Z is a matrix fo simulated vectors
    # y is the jth column of observed nominal responses
    yrep <- rep(0, dim(Z)[1])
    yrep[apply(Z, 1, max) < 0] <- 1
    yrep[yrep!=1] <- apply(Z[yrep!=1, ], 1, which.max) + 1
    
    probs <- as.vector(table(yrep)/dim(Z)[1])
    
    if(length(probs) < (dim(Z)[2]+1)) print("ERROR:No observations of one or more levels")
    
    Ez_nom <- matrix(NA, dim(Z)[2], dim(Z)[2]+1) 
    Ezzt_nom <- matrix(NA, dim(Z)[2], dim(Z)[2]+1)
    for(k in 1:(dim(Z)[2]+1)){
      Ez_nom[, k] <- apply(matrix(Z[yrep==k, ], nrow=sum(yrep==k)), 2, mean)
      Ezzt_nom[, k] <- apply(matrix(Z[yrep==k, ]^2, nrow=sum(yrep==k)), 2, mean)
    }
    list(probs, Ez_nom, Ezzt_nom)
  }
