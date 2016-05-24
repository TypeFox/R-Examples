direct_numerical_maximization <-
function(x, m,  delta, gamma, distribution_class, distribution_theta, DNM_limit_accuracy = 0.001, DNM_max_iter = 50, DNM_print = 2)
{
	
################################################################################################################################################################################################################################# Needed variables and functions ##############################################################################################
################################################################################################################################################################################  		
  discr_logL <- FALSE     ###### This algorithm hasn't been extendend for the use with the discrete likelihood, yet ############################################################	       
  discr_logL_eps <- 0.5   ###### This algorithm hasn't been extendend for the use with the discrete likelihood, yet ############################################################
  
  
  if (distribution_class == "pois" | distribution_class == "geom")
  {
    k <- 1
  }	
  if (distribution_class == "norm" | distribution_class == "genpois" | distribution_class == "bivariate_pois")
  {
    k <- 2
  }	
  
  
  DNM_log_n2w <- function(np)
  {
    wp <- log(np)
    return(wp)
  }
  
  
  DNM_exp_w2n <- function(wp)
  {
    np <- exp(wp)
    return(np)
  }
  
  
  DNM_logit_n2w <- function(np)
  {
    wp=log(np / (1 - np))
    return(wp)
  }
  
  
  DNM_invlogit_w2n <- function(wp)
  {
    np=exp(wp) / (1 + exp(wp))
    return(np)
  } 
  
  
  DNM_n2w <- function(m, gamma, distribution_class, distribution_theta)
  {	
    
    if (distribution_class == "pois")
    {	
      
      tlambda <- distribution_theta$lambda
      tlambda <- DNM_log_n2w(tlambda)
      tgamma <- NULL
      if (m > 1)
      {
        foo <- log(gamma / diag(gamma))
        tgamma <- as.vector(foo[!diag(m)])
      }
      parameter_vector <- c(tlambda,tgamma)
    }
    
    if (distribution_class == "norm")
    {	
      tmean <- distribution_theta$mean
      tsd <- distribution_theta$sd
      tsd <- DNM_log_n2w(tsd)
      tgamma <- NULL
      if (m > 1)
      {
        foo <- log(gamma / diag(gamma))
        tgamma <- as.vector(foo[!diag(m)])
      }
      parameter_vector <- c(tmean, tsd, tgamma)
    }
    
    if (distribution_class == "genpois")
    {	
      tlambda1 <- distribution_theta$lambda1
      tlambda1 <- DNM_log_n2w(tlambda1)
      tlambda2 <- distribution_theta$lambda2
      tlambda2 <- DNM_logit_n2w(tlambda2)
      tgamma <- NULL
      if (m > 1)
      {
        foo <- log(gamma / diag(gamma))
        tgamma <- as.vector(foo[!diag(m)])
      }
      parameter_vector <- c(tlambda1, tlambda2, tgamma)
    }
    return(parameter_vector)	
  }
   
   
  DNM_w2n <- function(parameter_vector, gamma, distribution_class, distribution_theta, m)
  {
    if (distribution_class == "pois")
    {
      distribution_theta$lambda <- DNM_exp_w2n(parameter_vector[1:m])
      gamma <- diag(m)
      epar <- exp(parameter_vector)
      if (m > 1)
      {
        gamma[!gamma] <- epar[((1 * m) + 1):((1 * m)+(m * m - m))]
        gamma <- gamma/apply(gamma,1,sum)
      }
      delta <- solve(t(diag(m)-gamma+1),rep(1,m))
    }
    
    if (distribution_class == "norm")
    {
      distribution_theta$sd <- DNM_exp_w2n(parameter_vector[(m + 1):(m + m)])
      distribution_theta$mean <- parameter_vector[1:m]
      gamma <- diag(m)
      epar <- exp(parameter_vector)
      if (m > 1)
      {
        gamma[!gamma] <- epar[((2 * m) + 1):(2 * m+(m * m - m))]
        gamma <- gamma/apply(gamma,1,sum)
      }
      delta <- solve(t(diag(m)-gamma+1),rep(1,m))
    }
    
    if (distribution_class == "genpois")
    {
      distribution_theta$lambda2 <- DNM_invlogit_w2n(parameter_vector[(m + 1):(m + m)])
      distribution_theta$lambda1 <- DNM_exp_w2n(parameter_vector[1:m])
      gamma <- diag(m)
      epar <- exp(parameter_vector)
      if (m > 1)
      {
        gamma[!gamma] <- epar[((2 * m)+1):(2 * m+(m * m - m))]
        gamma <- gamma / apply(gamma,1,sum)
      }
      delta <- solve(t(diag(m) - gamma + 1), rep(1,m))
    }
    return(list(delta=delta,
                gamma=gamma,
                distribution_theta = distribution_theta))
  }

  
  function_neglogL <- function(x, p, distribution_class, distribution_theta, m)
  {
    underflow_error_1=FALSE
    underflow_error_2=FALSE
    pos_logL_error=FALSE
    worse_logL_error=FALSE
    np <- DNM_w2n(parameter_vector = p, m = m, distribution_class = distribution_class, distribution_theta = distribution_theta)
    fb <-  try( forward_backward_algorithm(x = x, gamma = np$gamma, delta = np$delta, distribution_class = distribution_class, distribution_theta = np$distribution_theta,  discr_logL = discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)
    if (inherits(fb, "try-error"))
    {  
      underflow_error_1 <- TRUE
      fb$logL <- -Inf
    }
    neglogL <- (-1) * fb$logL
    
  return(neglogL)	
  }
  
################################################################################################################################################################################################################################# The algorithm ###############################################################################################################
################################################################################################################################################################################  
  
  if (distribution_class == "pois")
  {	
    if (any(distribution_theta$lambda == 0))
    { 
      distribution_theta$lambda <- ifelse(!distribution_theta$lambda == 0, distribution_theta$lambda, 1.1)
    }
    
    parameter_vector0 <- DNM_n2w(m = m, gamma = gamma, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    z <- nlm(f = function_neglogL, p = parameter_vector0, x = x, m = m, distribution_class = distribution_class,  distribution_theta = distribution_theta, print.level = DNM_print, gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
    logL <- - z$minimum
    
    par <- DNM_w2n(parameter_vector = z$estimate, m = m, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    delta <- par$delta
    
    gamma <- par$gamma
    
    distribution_theta <- par$distribution_theta
    
    estimated_mean_values <- distribution_theta$lambda
  }
  
  if (distribution_class == "norm")
  {	
    parameter_vector0 <- DNM_n2w(m = m, gamma = gamma, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    z <- nlm(f = function_neglogL, p = parameter_vector0, x = x, m = m, distribution_class = distribution_class, distribution_theta = distribution_theta, print.level = DNM_print, gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
    
    logL <- - z$minimum
    
    par <- DNM_w2n(parameter_vector = z$estimate, m = m, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    delta=par$delta
    
    gamma=par$gamma
    
    distribution_theta=par$distribution_theta
    
    estimated_mean_values <- distribution_theta$mean
  }
  
  if (distribution_class == "genpois")
  {	
    if (any(distribution_theta$lambda1 == 0))
    { 
      distribution_theta$lambda1 <- ifelse(!distribution_theta$lambda1 == 0, distribution_theta$lambda1, 0.1)
    }
    
    parameter_vector0 <- DNM_n2w(m = m, gamma = gamma, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    z <- nlm(f = function_neglogL, p = parameter_vector0, x = x, m = m, distribution_class = distribution_class,  distribution_theta = distribution_theta, print.level = DNM_print, gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
    
    logL <- - z$minimum
    
    par <- DNM_w2n(parameter_vector = z$estimate, m=m, distribution_class = distribution_class,distribution_theta = distribution_theta)
    
    delta <- par$delta
    gamma <- par$gamma
    distribution_theta <- par$distribution_theta
    estimated_mean_values <- distribution_theta$lambda1 / (1 - distribution_theta$lambda2)
  }
  
############################################################################################################################################################################################################################### Accessing AIC and BIC for the trained HMM ####################################################################################
################################################################################################################################################################################
  
  AIC <- AIC_HMM(logL = logL, m = m, k = k) 
  BIC <- BIC_HMM(size = length(x), logL = logL, m = m, k = k) 
  
############################################################################################################################################################################################################################### Return results ################################################################################################################
################################################################################################################################################################################  
  return(list(x = x, 
              m = m, 
              logL = logL, 
              AIC = AIC, 
              BIC = BIC, 
              delta = delta, 
              gamma = gamma, 
              distribution_class = distribution_class, 
              distribution_theta = distribution_theta, 
              estimated_mean_values=estimated_mean_values))
}
