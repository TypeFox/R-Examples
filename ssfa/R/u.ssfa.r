# Jodrow et al (1982)    
u.ssfa <- function(object, ...) {
  if(object$rho==0)
  {
    sigmau2   <- object$sigmau2
    sigmav2   <- object$sigmav2
    fun <- object$fun
    sc <- object$sc
    
    if (fun == "hnormal") {
      mu_i <- -sc * residuals.ssfa(object) * sigmau2 / (sigmau2 + sigmav2)
      sigma_i <- sqrt(sigmau2 * sigmav2 / (sigmau2 + sigmav2))
      
    }
    u <- mu_i + sigma_i*(dnorm(-mu_i/sigma_i)/pnorm(mu_i/sigma_i))
  }
    
    
  if(object$rho!=0)
  {     
  if(object$sigmau2_dmu<=0)
  {
    u <- matrix(0,ncol=1, nrow=length(object$y))
  }
  else
  {  
  sigmau2   <- object$sigmau2_dmu
  sigmav2   <- object$sigmav2
  fun <- object$fun
  sc <- object$sc
    
  if (fun == "hnormal") {
    mu_i <- -sc * residuals.ssfa(object) * sigmau2 / (sigmau2 + sigmav2)
    sigma_i <- sqrt(sigmau2 * sigmav2 / (sigmau2 + sigmav2))
  
  }
  u <- mu_i + sigma_i*(dnorm(-mu_i/sigma_i)/pnorm(mu_i/sigma_i))
  }
  }
  return(u)
}