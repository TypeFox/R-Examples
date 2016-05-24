"epi.empbayes" <- function(obs, pop){
	# gamma: mean of rate
  # phi: variance of rate
        
  gamma <- (sum(obs)) / (sum(pop))
  rate <- obs / pop
  sum.pop <- sum(pop)
        
  phi.left <- sum(pop * (rate - gamma)^2) / sum.pop
  phi.right <- gamma / mean(pop)
  phi <- phi.left - phi.right
        
  # The convention is that phi = 0 whenever the above expression is negative.
  phi <- ifelse(phi < 0, 0, phi)
  emp <- ((phi * (rate - gamma)) / (phi + (gamma / pop))) + gamma
        
  # gamma = nu / alpha
  # phi = nu / alpha^2
        
  alpha <- gamma / phi
  nu <- gamma^2 / phi
  inv.nu <- 1 / nu
                
  rval <- data.frame(gamma, phi, alpha, nu, inv.nu)
  names(rval) <- c("gamma (mean)", "phi (variance)", "alpha (shape)", "nu (scale)", "inv.nu (rate)")
  unlist(rval)
}
    
    
