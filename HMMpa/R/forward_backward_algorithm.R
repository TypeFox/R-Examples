forward_backward_algorithm <-
function(x, delta, gamma, distribution_class, distribution_theta, discr_logL = FALSE, discr_logL_eps = 0.5)
{

################################################################################################################################################################################################################################# Needed variables and functions ##############################################################################################
################################################################################################################################################################################  	
  size <- length(x)
  m <- length(delta)
  fb_svtpoue <- fb_small_value_to_prevent_overflow_and_underflow_errors <- 4.940656e-142   #6.716723e-269   

  function_discr_log_L_p_norm <- function(x, mean, sd, discr_logL_eps)
  {
    foo <- pnorm((x + discr_logL_eps), mean = mean, sd = sd) - pnorm((x - discr_logL_eps), mean = mean, sd = sd)	
  return(foo)
  }    
    
################################################################################################################################################################################################################################# Preparation for calcualtion of log_alphas and log_betas    ##################################################################
################################################# (Calcualtion of probabilities)    ############################################################################################
################################################################################################################################################################################  

  if (distribution_class == "pois")
  {	
	#### old
    #probabilities <- matrix(0, ncol = m, nrow = size)
    #for (i in 1:size)
    #{
    #  for (j in 1:m)
    #  {
    #    probabilities[i,j] <- dpois(x[i], lambda = distribution_theta$lambda[j])
    #  } 
    #}	
  
    ### old_2
    # 
  	# probabilities <-  matrix(x, ncol = m, nrow = size)
    # probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN = dpois, lambda = distribution_theta$lambda))
  
    probabilities <- outer(X = x, Y = distribution_theta$lambda, FUN = dpois)
  }
  
  
  if (distribution_class == "geom")
  {  
    probabilities <-  matrix(x, ncol = m, nrow = size)
    probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN = dgeom, prob = distribution_theta$prob))
  }
  

  if (distribution_class == "genpois")
  {
    probabilities <-  matrix(x, ncol = m, nrow = size)
    probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN = dgenpois, lambda1 = distribution_theta$lambda1, lambda2 = distribution_theta$lambda2))
  }
  
  
  if (distribution_class == "norm" & discr_logL == FALSE)
  {
    probabilities <-  matrix(x, ncol = m, nrow = size)
    probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN = dnorm, mean = distribution_theta$mean, sd = distribution_theta$sd))
  }
  
  
  if (distribution_class == "norm" & discr_logL == TRUE)
  {  
    probabilities <-  matrix(x, ncol = m, nrow = size)
    probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN =  function_discr_log_L_p_norm, mean = distribution_theta$mean, sd = distribution_theta$sd, discr_logL_eps = discr_logL_eps))
  }
  
      
  if (distribution_class == "bivariate_pois")
  {   
    size <- length(x[,1])
    
    probabilities <- matrix(0, ncol = m, nrow = size)
    for (i in 1:size)
    {
      for (j in 1:m)
      {
        probabilities[i,j] <- dpois(x[i,1], lambda=distribution_theta$lambda_1[j]) * dpois(x[i,2],lambda=distribution_theta$lambda_2[j])
      } 
    }	
  }
   ################################################################################################################################################################################################################################ Negative probabilitiesabilities etc. must be caused by numerical over- and or underflow #################################### ################################################# -> setting those results to a very small number (fb_svtpoue) (Caution: small manipulation) ##################################
###############################################################################################################################################################################
  
    probabilities <- ifelse(!is.na(probabilities), probabilities, fb_svtpoue)
    probabilities <- ifelse(!probabilities <= 0, probabilities, fb_svtpoue) 	
    probabilities <- ifelse(!probabilities > 1, probabilities, fb_svtpoue) 	
    probabilities <- ifelse(!probabilities == Inf, probabilities, fb_svtpoue)	
    probabilities <- ifelse(!probabilities == -Inf, probabilities, fb_svtpoue) 
    
################################################################################################################################################################################################################################# Calculating log_alphas         ##############################################################################################
################################################################################################################################################################################  
  
  log_alpha <- matrix(fb_svtpoue, nrow = size, ncol = m)
  foo <- delta * probabilities[1,]
  sum_of_foo <- sum(foo) + fb_svtpoue
  scaled_logL <- log(sum_of_foo)
  foo <- foo / sum_of_foo
  log_alpha[1,] <- scaled_logL + log(foo)
  for (i in 2:size)
  {
    foo <- foo %*% gamma * probabilities[i,]
    sum_of_foo <- sum(foo) + fb_svtpoue
    scaled_logL <- scaled_logL + log(sum_of_foo)
    foo <- foo / sum_of_foo
    log_alpha[i,] <- scaled_logL + log(foo)
  }
  ################################################# Calculating log_L via alpha_T ###########################################################################################
  logL_calculated_with_alpha_T <- scaled_logL
  
  
################################################################################################################################################################################################################################# Calculating log_betas         ###############################################################################################
################################################################################################################################################################################
   
  
  log_beta <- matrix(fb_svtpoue, nrow = size, ncol = m)
  log_beta[size,] <- rep(0,m)
  foo <- rep(1 / m, m)
  scaled_logL <- log(m)
  for (i in (size-1):1)
  {
    foo <- gamma %*% (probabilities[i+1,] * foo)
    log_beta[i,] <- log(foo) + scaled_logL
    sum_of_foo <- sum(foo) + fb_svtpoue
    foo <- foo / sum_of_foo
    scaled_logL <- scaled_logL + log(sum_of_foo)
  }
  ################################################# Calculating log_L via beta_1 ###############################################################################################
  logL_calculated_with_beta_1 <- scaled_logL
  
  
################################################################################################################################################################################################################################# Calculating log_L via alpha_middle_t and beta_middle_t         ##############################################################
################################################################################################################################################################################
     
  logL_calculated_with_alpha_t_and_beta_t <- 0
  middle_t <- round(size / 2)
  c <- max(log_alpha[middle_t,])
  for (i in 1:m)
  { 
  	logL_calculated_with_alpha_t_and_beta_t <- logL_calculated_with_alpha_t_and_beta_t + exp( log_alpha[middle_t,i] + log_beta[middle_t,i] - c)
  }
  logL_calculated_with_alpha_t_and_beta_t <- log(logL_calculated_with_alpha_t_and_beta_t) + c
  
  logL <- logL_calculated_with_alpha_t_and_beta_t
  logL_calculation <- "logL calculated with alpha_middle_t and beta_middle_t"
  if (logL == -Inf | logL == Inf | is.na(logL))
  {   
    logL_calculation <- "logL calculated with alpha_T"
    logL <- logL_calculated_with_alpha_T	
  }
  if (logL == -Inf | logL == Inf | is.na(logL))
  {
    logL <- logL_calculated_with_beta_1
    logL_calculation <- "logL calculated with beta_1"
  }

############################################################################################################################################################################################################################### Return results ################################################################################################################
################################################################################################################################################################################  
return(list(log_alpha = log_alpha, 
            log_beta = log_beta, 
            logL = logL, 
            logL_calculation = logL_calculation))			
}
