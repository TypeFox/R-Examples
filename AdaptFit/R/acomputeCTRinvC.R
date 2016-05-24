########## R function acompute.CTRinvC #############

# Computes the inverse residual covariance matrix of model coefficients
# for linear mixed models


"acompute.CTRinvC" <-
  function (X, Z, RR.inv, individual, rho = 0) 
{
  C <- cbind(X, Z)
  CTRinvC <- 0 * t(C) %*% C
  if (length(individual) > 1) {
    end.indiv <- (1:length(individual))[diff(individual) != 
                                        0]
    strt.indiv <- c(1, end.indiv + 1)
    end.indiv <- c(end.indiv, length(individual))
    n.indiv <- length(end.indiv)
    indiv.inds <- NULL
    indiv.times <- NULL
    for (i in 1:n.indiv) {
      indiv.inds <- (strt.indiv[i]:end.indiv[i])
      indiv.times <- (1:length(indiv.inds))
      R.inv <- solve(rho^abs(outer(indiv.times, indiv.times, 
                                   "-")))
      if (ncol(R.inv) > 1) {
        CTRinvC <- CTRinvC + t(C[indiv.inds, ]) %*% R.inv %*% 
          C[indiv.inds, ]
      }
      else {
        CTRinvC <- CTRinvC + as.numeric(R.inv) * outer(C[indiv.inds, 
                                                         ], C[indiv.inds, ])
      }
    }
  }
  else 
      if(!is.null(RR.inv))
         CTRinvC <- t(C) %*% RR.inv%*%C
  else
     CTRinvC <- t(C)%*%C
    
    
  return(CTRinvC)
}
