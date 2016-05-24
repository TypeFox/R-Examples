HMM_decoding <-
function(x, m, delta, gamma, distribution_class, distribution_theta, decoding_method = "global", discr_logL = FALSE, discr_logL_eps = 0.5)
{

################################################################################################################################################################################################################################# Needed variables and functions ##############################################################################################
################################################################################################################################################################################  	
  svtpoue <- small_value_to_prevent_overflow_and_underflow_errors <- 4.940656e-142
  
  size <- length(x)
      
################################################################################################################################################################################################################################# Preparation for calcualtion of local most probable states ##################################################################
################################################################################################################################################################################  
  
  if (distribution_class == "pois")
  {
    distribution_means <- distribution_theta$lambda		
  }
  
  if (distribution_class == "geom")
  {
    distribution_means <- 1 / distribution_theta$prob$prob					
  }
  
  if (distribution_class == "genpois")
  {	
    distribution_means <- distribution_theta$lambda1 / (1 - distribution_theta$lambda2) 
  }
  
  if(distribution_class == "norm" & discr_logL == FALSE)
  {	
    distribution_means <- distribution_theta$mean
  }
  
  if(distribution_class == "norm" & discr_logL == TRUE)
  {	
    distribution_means <- distribution_theta$mean
  }

################################################################################################################################################################################################################################# Calculation of the most probable states  ####################################################################################
################################################################################################################################################################################  
  
 if (decoding_method == "global") 
 {   
 	decoding <- Viterbi_algorithm(x = x, m = m, delta = delta, gamma = gamma, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps)$decoding
 }
 
 if (decoding_method == "local") 
 {   
 	decoding <- local_decoding_algorithm(x = x, m = m, delta = delta, gamma = gamma, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL = discr_logL, discr_logL_eps = discr_logL_eps)$decoding
 }
 
################################################################################################################################################################################################################################# Connection of local most probable state with the mean of the state-dependend distributions   ###############################
################################################################################################################################################################################  
  
 decoding_distr_means <- distribution_means[decoding]
 
############################################################################################################################################################################################################################### Return results ################################################################################################################
################################################################################################################################################################################
  
return(list(decoding_method = decoding_method, 
			      decoding = decoding, 
			      decoding_distr_means = decoding_distr_means))
}
