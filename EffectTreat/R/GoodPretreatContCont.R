GoodPretreatContCont <- function(T0T0, T1T1, Delta, T0T1=seq(from=0, to=1, by=.01)){

sigma_delta_T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * T0T1)
rho2_min <- 1 - (Delta/sigma_delta_T)   
T0T1 <- T0T1[rho2_min >=-1 & rho2_min <=1]
sigma_delta_T <- sigma_delta_T[rho2_min >=-1 & rho2_min <=1] 
rho2_min <- rho2_min[rho2_min >=-1 & rho2_min <=1] 

fit <- 
  list(T0T1=T0T1, Sigma.Delta.T=sigma_delta_T, Rho2.Min = rho2_min, Call=match.call())   

class(fit) <- "GoodPretreatContCont"
fit
}
