# - gamma_v and its derivative for the dynamic model 

Gamma_v_f <- function(Rho, T) {
  rho_power <- Rho^(0:(T-1))
  rho_d <- c(0,(Rho^(0:(T-2)))*(1:(T-1)) )
  
  list(Gamma_v= rho_power %*% t(rho_power),
       d.Gamma_v = rho_power %*% t(rho_d) +
                  rho_d %*% t(rho_power))
 }
