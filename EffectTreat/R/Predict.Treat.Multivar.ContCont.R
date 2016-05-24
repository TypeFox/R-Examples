Predict.Treat.Multivar.ContCont <- function(Sigma_TT, Sigma_TS, Sigma_SS, Beta, S, mu_S, 
                                            T0T1=seq(-1, 1, by=.01)) { 

  Results <- T0T1_results <- var_Delta.T_given_S_results <- NULL
  
  for (i in 1: length(T0T1)) {   
    
    Sigma_TT_hier <- Sigma_TT
    Sigma_TT_hier[1,2] <- Sigma_TT_hier[2, 1] <- (sqrt(Sigma_TT[1,1])*sqrt(Sigma_TT[2,2])) * T0T1[i]
    
    a_1 <- matrix(c(-1, 1), nrow=1)
    
    Sigma_temp_a <- 
      rbind(cbind(Sigma_TT_hier, Sigma_TS), 
            cbind(t(Sigma_TS), Sigma_SS))
    Cor_c <- cov2cor(Sigma_temp_a)
    a_1 <- matrix(c(-1, 1), nrow=1)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)  
    
    
    if (Min.Eigen.Cor > 0) {
      
      Pred_T <- 
        Beta + a_1 %*% Sigma_TS %*% solve(Sigma_SS) %*% (S - mu_S)
      
      Results <- rbind(Results, Pred_T)
      T0T1_results <- rbind(T0T1_results, T0T1[i])
      var_Delta.T_given_S <- 
        (a_1 %*% Sigma_TT_hier %*% t(a_1)) - 
        (a_1 %*% Sigma_TS %*% solve(Sigma_SS) %*% t(Sigma_TS) %*% t(a_1)) 
      var_Delta.T_given_S_results <- rbind(var_Delta.T_given_S_results, var_Delta.T_given_S)
    }
  }
  
fit <- 
    list(Pred_T=Results[1], Var_Delta.T_S=as.numeric(t(var_Delta.T_given_S_results)), T0T1=as.numeric(t(T0T1_results)), Call=match.call())
  
  class(fit) <- "Predict.Treat.Multivar.ContCont"
  fit
}