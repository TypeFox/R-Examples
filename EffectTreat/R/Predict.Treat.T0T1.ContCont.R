Predict.Treat.T0T1.ContCont <- function(x, S, Beta, SS, mu_S, T0T1, alpha=.05) { 
  
  if (T0T1 < min(x$Pos.Def$T0T1) | T0T1 > max(x$Pos.Def$T0T1)){
    stop("The specified value T0T1 should be within the range of values for which valid solutions were obtained
         (examine x$Pos.Def$T0T1, where x is the fitted object of class PCA.ContCont)")
  } 
  
  user_req_Var_Delta.T <- x$GoodSurr$Sigma.Delta.T[which.min(abs(x$Pos.Def$T0T1-T0T1))]
  PCA_val <- x$PCA[which.min(abs(x$Pos.Def$T0T1-T0T1)) ]
  
  g_S_j <- Beta + (sqrt(user_req_Var_Delta.T / SS)) * PCA_val * (S - mu_S)
  
  var_Delta.T <-
    user_req_Var_Delta.T 
  
  var_Delta.T_given_S <- 
    user_req_Var_Delta.T * (1-(PCA_val**2))
  
  crit_val <- qnorm(c(alpha/2), mean=0, sd=1, lower.tail=FALSE)
  
  CIs_LB <- (-crit_val*sqrt(var_Delta.T_given_S))+g_S_j
  CIs_UB <- (crit_val*sqrt(var_Delta.T_given_S))+g_S_j
  
  fit <- 
    list(Pred_T=g_S_j[1], Var_Delta.T=var_Delta.T, T0T1=T0T1, CI_low=CIs_UB, CI_high=CIs_LB, Var_Delta.T_S=var_Delta.T_given_S, alpha=alpha,
         Call=match.call())
  
  class(fit) <- "Predict.Treat.T0T1.ContCont"
  fit
  
  }  
