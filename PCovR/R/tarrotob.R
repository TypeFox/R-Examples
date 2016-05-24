tarrotob <-
function(F1,W){
  conv <- 1e-6
  m <- nrow(F1)
  r <- ncol(F1)
  TA  <- diag(r)
  PHI <- diag(r)
  A   <- F1 %*% TA
  
  A2 <- (A-W*A)^2
  f1 <- sum(A2)
  
  f1old <-  f1
  
  iter <-  1
  again <-  1
  while (again ==1){
    for (g in 1:r){
      for (h in 1:r){
        if (g!=h){
          agg <-0
          bgg <-0
          bgh <-0
          for (j in 1:m){
            agg <- agg + ((1-W[j,g]) * A[j,g]^2)
            bgg <- bgg + ((1-W[j,h]) * A[j,g]^2)
            bgh <- bgh + ((1-W[j,h]) * A[j,g] * A[j,h])
          }
          
          delta <- (bgh-PHI[g,h]*agg)/(agg+bgg)
          gamma <- sqrt(delta^2 + 2*PHI[g,h]*delta+1)
          
          for (j in 1:m){
            A[j,h] <- (-1.0 * delta * A[j,g])+ A[j,h]
            A[j,g] <- gamma * A[j,g]
          }
          
          for (j in 1:r){
            TA[j,h] <- (-1.0 * delta * TA[j,g])+ TA[j,h]
            TA[j,g] <- gamma * TA[j,g]
          }
          
          for (j in 1:r){
            PHI[g,j] <- (PHI[g,j]+ delta*PHI[h,j])/gamma
          }
          PHI[g,g] <- 1
          for (j in 1:r){
            PHI[j,g] = PHI[g,j]
          }
        }
        
      }  
    }
    A2 <- (A-W*A)^2
    f1 <- sum(A2)
    
    #stop criteria
    if (f1 < .1*conv || abs (f1 - f1old) < conv || iter > 50){
      again <- 0
    }
    iter <- iter + 1
  }
  ST <- A %*% PHI
  T <- TA
  res <- list(T=T,A=A)
  return(res)
}
