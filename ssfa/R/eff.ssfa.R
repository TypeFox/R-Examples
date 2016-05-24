eff.ssfa <- function(object, ...) {
  
  if(object$rho==0)
  {
    sigmau2   <- object$sigmau2
    sigmav2   <- object$sigmav2
    sc <- object$sc
    
    if (sc==1) #production
    {
      mu_i <-  -sc * residuals.ssfa(object) * sigmau2 / (sigmau2 + sigmav2)
      sigma_i <- sqrt(sigmau2) * sqrt(sigmav2) / sqrt(sigmau2 + sigmav2)
      eff <- (1 - pnorm(sc * sigma_i - mu_i / sigma_i)) / (1 - pnorm(-mu_i / sigma_i))* exp(-sc * mu_i + 0.5 * sigma_i^2)
    }
    else
    {
      mu_i <-  -sc * residuals.ssfa(object) * sigmau2 / (sigmau2 + sigmav2)
      sigma_i <- sqrt(sigmau2) * sqrt(sigmav2) / sqrt(sigmau2 + sigmav2)
      eff <- 1/((1 - pnorm(sc * sigma_i - mu_i / sigma_i)) / (1 - pnorm(-mu_i / sigma_i))* exp(-sc * mu_i + 0.5 * sigma_i^2))
    }
  }
  
  if(object$rho!=0)
  { 
  if(object$sigmau2_dmu<=0)
  {
    eff <- matrix(1,ncol=1, nrow=length(object$y))
  }
  else
  {
  sigmau2   <- object$sigmau2_dmu
  sigmav2   <- object$sigmav2
  sc <- object$sc
  
  if (sc==1) #production
    {
    mu_i <-  -sc * residuals.ssfa(object) * sigmau2 / (sigmau2 + sigmav2)
    sigma_i <- sqrt(sigmau2) * sqrt(sigmav2) / sqrt(sigmau2 + sigmav2)
    eff <- (1 - pnorm(sc * sigma_i - mu_i / sigma_i)) / (1 - pnorm(-mu_i / sigma_i))* exp(-sc * mu_i + 0.5 * sigma_i^2)
    }
  else
  {
    mu_i <-  -sc * residuals.ssfa(object) * sigmau2 / (sigmau2 + sigmav2)
    sigma_i <- sqrt(sigmau2) * sqrt(sigmav2) / sqrt(sigmau2 + sigmav2)
    eff <- 1/((1 - pnorm(sc * sigma_i - mu_i / sigma_i)) / (1 - pnorm(-mu_i / sigma_i))* exp(-sc * mu_i + 0.5 * sigma_i^2))
  }
  }
  }
  return(eff)
  }





