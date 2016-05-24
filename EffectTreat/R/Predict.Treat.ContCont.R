Predict.Treat.ContCont <- function(x, S, Beta, SS, mu_S) { 
  
Sigma.Delta.T <- x$GoodSurr$Sigma.Delta.T
PCA_vec <- x$PCA
T0T1 <- x$Pos.Def$T0T1  

g_S_j <- 
  Beta + (sqrt(Sigma.Delta.T / SS)) * PCA_vec * (S - mu_S)
 
var_Delta.T <-
  Sigma.Delta.T 

var_Delta.T_given_S <- 
      Sigma.Delta.T * (1-(PCA_vec**2))

fit <- 
    list(Pred_T=g_S_j[1], Var_Delta.T=var_Delta.T, T0T1=T0T1, PCA=PCA_vec, Var_Delta.T_S=var_Delta.T_given_S, Call=match.call())
  
  class(fit) <- "Predict.Treat.ContCont"
  fit
  
}