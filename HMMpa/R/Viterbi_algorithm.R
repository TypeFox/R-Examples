Viterbi_algorithm <-
function(x, m, delta, gamma, distribution_class, distribution_theta, discr_logL = FALSE, discr_logL_eps = 0.5)
{

################################################################################################################################################################################################################################# Needed variables and functions ##############################################################################################
################################################################################################################################################################################  	
  svtpoue <- small_value_to_prevent_overflow_and_underflow_errors <- 4.940656e-142
  
  size <- length(x)
  
  function_discr_log_L_p_norm <- function(x, mean, sd, discr_logL_eps)
  {
    foo <- pnorm((x + discr_logL_eps), mean = mean, sd = sd) - pnorm((x - discr_logL_eps), mean = mean, sd = sd)	
  return(foo)
  }    
    
################################################################################################################################################################################################################################# Preparation for calcualtion of global most probable states ##################################################################
################################################# (Calcualtion of probabilities)    ############################################################################################
################################################################################################################################################################################  
  
  if (distribution_class == "pois") 
  {
    distribution_means <- distribution_theta$lambda		
    
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
  	# probabilities <- matrix(x, ncol = m, nrow = size)
    # probabilities <- t(apply(X = probabilities, MARGIN = 1, FUN = dpois, lambda = distribution_theta$lambda))
    
     probabilities <- outer(X = x, Y = distribution_theta$lambda, FUN = dpois)
  }
  
  if (distribution_class == "geom") 
  {
    distribution_means <- 1 / distribution_theta$prob$prob				
	probabilities <- matrix(x, ncol = m, nrow = size)
    probabilities <- t(apply(X = probabilities, MARGIN = 1, FUN = dgeom, prob = distribution_theta$prob))	
  }
  
  if (distribution_class == "genpois") 
  {	 
    distribution_means <- distribution_theta$lambda1 / (1 - distribution_theta$lambda2) 
    probabilities <- matrix(x, ncol = m, nrow = size)
    probabilities <- t(apply(X = probabilities, MARGIN = 1, FUN = dgenpois, lambda1 = distribution_theta$lambda1, lambda2 = distribution_theta$lambda2))
  }
  
  if (distribution_class == "norm" & discr_logL == FALSE)
  {	
    distribution_means <- distribution_theta$mean
    probabilities <- matrix(x, ncol = m, nrow = size)
    probabilities <- t(apply(X = probabilities, MARGIN = 1, FUN = dnorm, mean = distribution_theta$mean, sd = distribution_theta$sd))
  }
  
  if (distribution_class == "norm" & discr_logL == TRUE)
  {	
    distribution_means <- distribution_theta$mean
    probabilities <- matrix(x, ncol = m, nrow = size)
    probabilities <- t(apply(X = probabilities, MARGIN = 1, FUN = function_discr_log_L_p_norm, mean = distribution_theta$mean, sd = distribution_theta$sd, discr_logL_eps = discr_logL_eps))
  }

################################################################################################################################################################################################################################ Negative probabilitiesabilities etc. must be caused by numerical over- and or underflow ################################################ -> setting those results to a very small number (fb_svtpoue) (Caution: small manipulation) ##################################
###############################################################################################################################################################################

    probabilities <- ifelse(!is.na(probabilities), probabilities, svtpoue)
    probabilities <- ifelse(!probabilities <= 0, probabilities, svtpoue) 		
    probabilities <- ifelse(!probabilities == Inf, probabilities, svtpoue)	
    probabilities <- ifelse(!probabilities == -Inf, probabilities, svtpoue) 

################################################################################################################################################################################################################################# The algorithm     ###########################################################################################################
################################################################################################################################################################################  
  
  omega <- matrix(svtpoue, ncol = m, nrow = size)
  foo <- delta * probabilities[1,]
  omega[1,] <- foo / sum(foo)
  for(i in 2:size)
  {
    foo <- apply(omega[i-1,] * gamma, 2, max) * probabilities[i,]
    omega[i,] <- foo / sum(foo)
  }
  
  decoding <- numeric(size)
  decoding[size] <- which.max(omega[size,])
  for(i in (size - 1):1)
  {
    decoding[i] <- which.max(gamma[,decoding[i+1]] * omega[i,])
  }


############################################################################################################################################################################################################################### Return results ################################################################################################################
################################################################################################################################################################################
  
return(list(omega=omega, 
            decoding=decoding))
}
