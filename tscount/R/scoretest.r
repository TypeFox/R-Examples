scoretest <- function(Score, G, G1, r, ...){
#Computes score test statistic 
##############################
  R <- seq(along=numeric(r))
  s <- length(Score)-r #number of parameters of the restricted model
  S <- seq(along=numeric(s))
  Score_2 <- Score[s+R]
  G_11 <- G[S, S]
  G_12 <- G[S, s+R]
  G_21 <- G[s+R, S]
  G_22 <- G[s+R, s+R]
  G1_11 <- G1[S, S]
  G1_12 <- G1[S, s+R]
  G1_21 <- G1[s+R, S]
  G1_22 <- G1[s+R, s+R]
  G_11_inv_temp <- invertinfo(G_11, ...)
  if(!is.null(G_11_inv_temp$error_message)) return(list(error_message=G_11_inv_temp$error_message))
  G_11_inv <- G_11_inv_temp$vcov
  Sigma <- G1_22 - G_21%*%G_11_inv%*%G1_12 - G1_21%*%G_11_inv%*%G_12 + G_21%*%G_11_inv%*%G1_11%*%G_11_inv%*%G_12
  Sigma_inv_temp <- invertinfo(Sigma, ...)  
  if(!is.null(Sigma_inv_temp$error_message)) return(list(error_message=Sigma_inv_temp$error_message))
  Sigma_inv <- Sigma_inv_temp$vcov
  result <- list(test_statistic=(t(Score_2) %*% Sigma_inv %*% Score_2)[1, 1])
  return(result)
}
