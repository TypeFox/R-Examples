calculateGaussianBc <-
function(model, sigma.estimated, analytic) {
  # A function that calculates the analytic representation of the bias 
  # corrections in linear mixed models, see Greven & Kneib (2010).
  #
  # Args: 
  #   model    = From getAllModelComponents()
  #   sigma.estimated = If sigma is estimated. This only is used for the 
  #                     analytical version of Gaussian responses.
  #   analytic = FALSE if the numeric hessian of the (restricted) marginal log-
  #              likelihood from the lmer optimization procedure should be used.
  #              Otherwise (default) TRUE, i.e. use a analytical version that 
  #              has to be computed.
  #
  # Returns:
  #   df = Bias correction (i.e. degrees of freedom) for a linear mixed model.
  #           
  C     <- model$C
  B     <- model$B
  e     <- model$e
  A     <- model$A
  tye   <- model$tye
  V0inv <- model$V0inv
   
  if(analytic) {
    for (j in 1:length(model$theta)) {
      Wj     <- model$Wlist[[j]]
      eWje   <- model$eWelist[[j]]
      C[j, ] <- as.vector((e %*% Wj) %*% A - eWje * e/(2 * tye))
      for (k in j:length(model$theta)) {
          Wk <- model$Wlist[[k]]
          eWke   <- model$eWelist[[k]]
          if (!model$isREML) {
            B[j, k] <- B[k, j] <-  - tye * 
              sum(t(Wk %*% V0inv) * (Wj %*% V0inv))/(2 * model$n) - 
              eWje * eWke/(2 * tye) + 
              as.numeric(e %*% Wk %*% (A %*% (Wj %*% e)))
          } else {
            B[j, k] <- B[k, j] <- - tye * 
              sum(t(Wk %*% A) * (Wj %*% A))/(2*(model$n - ncol(model$X))) - 
              eWje * eWke/(2 * tye) + 
              as.numeric(e %*% Wk %*% (model$A %*% (Wj %*% e)))
          }
      }
    }
  } else {
    if (!model$isREML) {
      np <- model$n
    } else {
      np <- model$n - ncol(model$X)
    }
    
    for (j in 1:length(model$theta)) {
        Wj     <- model$Wlist[[j]]
        eWje   <- model$eWelist[[j]]
        C[j, ] <- 2 * np / tye * as.vector((e %*% Wj) %*% A - eWje * e/tye)
    }
  }
  
  Rchol   <- chol(B)
  L1      <- backsolve(Rchol, C, transpose = TRUE)
  Lambday <- backsolve(Rchol, L1)

  df <- model$n - sum(diag(A))
  for (j in 1:length(model$theta)) {
      df <- df + sum(Lambday[j,] %*% (A %*% (model$Wlist[[j]] %*% e)))  
  }
  
  if (sigma.estimated) {
    df <- df + 1
  }
  
  return(df)
}
