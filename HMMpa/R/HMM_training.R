HMM_training <-
function(x, distribution_class, min_m = 2, max_m = 6, 
         n = 100, training_method = "EM", discr_logL = FALSE, discr_logL_eps = 0.5, 
         Mstep_numerical = FALSE, dynamical_selection = TRUE, 
         BW_max_iter = 50, BW_limit_accuracy = 0.001, BW_print = TRUE,
         DNM_max_iter = 50, DNM_limit_accuracy = 0.001, DNM_print = 2)
{

################################################################################################################################################################################################################################# Needed variables ############################################################################################################
################################################################################################################################################################################  

		
  par_var <- 0.5 
  how_many_HMMs <- length(min_m:max_m) + (min_m - 1)
  
  ############### define the arrays to save the estimates for each m-state-HMM    #############################################################################################
  list_of_all_initial_parameters <- pairlist(NA)
  list_of_all_trained_HMMs <- pairlist(NA)
  list_of_all_logLs_for_each_HMM_with_m_states <- rep(NA, times = min_m - 1)
  list_of_all_AICs_for_each_HMM_with_m_states <- rep(NA, times = min_m - 1)
  list_of_all_BICs_for_each_HMM_with_m_states <- rep(NA, times = min_m - 1)
  
     ############################################################################################################################################################################################################################### The training        #########################################################################################################
################################################################################################################################################################################
 
  for (i in min_m:max_m) 
  {	
  	############### print each step to see for which m the HMM is trained    ###############################################################################################
  	print(paste('HMM with m =', toString(i)))
    
    ############################################################################################################################################################################
    ############### intial parameter training/setting            ###############################################################################################################
    ############################################################################################################################################################################
    if (dynamical_selection == FALSE) 
    {
      temp_initial_parameters <- initial_parameter_training(n = n, x = x, m = i, distribution_class = distribution_class,  discr_logL= discr_logL, discr_logL_eps = discr_logL_eps)
    }			
    if (dynamical_selection == TRUE) 
    {
      if (i == min_m) 
      {
        temp_initial_parameters <- initial_parameter_training(n = n, x = x, m = i, distribution_class = distribution_class, discr_logL = discr_logL, discr_logL_eps = discr_logL_eps)
      }		
      if (i > min_m) 
      {	
        temp_initial_parameters$m <- i				
        temp_initial_parameters$delta <- c(temp_trained_HMM_regarding_distribution_theta$delta , 1 / i)
        temp_initial_parameters$delta <- temp_initial_parameters$delta / sum(temp_initial_parameters$delta) 				
        temp_initial_parameters$gamma <- matrix(0.2 / i, nrow = i, ncol = i)
        temp_initial_parameters$gamma[i,i] <- temp_initial_parameters$gamma[i,i] + 0.8
        for (j in 1:(i - 1)) 
        {
          temp_initial_parameters$gamma[j,1:(i - 1)] <- temp_trained_HMM_regarding_distribution_theta$gamma[j,]
        }
        for (j in 1:i) 
        {
          temp_initial_parameters$gamma[j,] <- temp_initial_parameters$gamma[j,] / sum(temp_initial_parameters$gamma[j,])
        }				
        temp_distance <- x
        for (j in 1:length(x)) 
        {
          temp_distance[j] <- min(abs(x[j] - temp_trained_HMM_regarding_distribution_theta$estimated_mean_values))
        }
        temp_next_estimated_mean_value <- x[which.max(temp_distance)]
        temp_initial_parameters$E <- c(temp_trained_HMM_regarding_distribution_theta$estimated_mean_values, temp_next_estimated_mean_value)				
        if (distribution_class == "pois") 
        {
          temp_initial_parameters$distribution_theta <- list(lambda=temp_initial_parameters$E)
        }
        if (distribution_class == "genpois") 
        {
          temp_initial_parameters$distribution_theta <- list(lambda1 = temp_initial_parameters$E * (1 - par_var) , lambda2 = rep(par_var, times = i))
        }
        if (distribution_class == "norm") 
        {
          temp_initial_parameters$distribution_theta <- list(mean = temp_initial_parameters$E, sd = rep(sd(x) / i, times = i))
        }
        if (distribution_class == "geom") 
        {
          temp_initial_parameters$distribution_theta <- list(prob = 1 / temp_initial_parameters$E)
        }				
        temp_initial_parameters$logL <- forward_backward_algorithm(x = x, delta = temp_initial_parameters$delta, gamma = temp_initial_parameters$gamma, distribution_class = distribution_class, distribution_theta = temp_initial_parameters$distribution_theta,  discr_logL = discr_logL, discr_logL_eps = discr_logL_eps)$logL
      }
    }		
    list_of_all_initial_parameters[[i]] <- temp_initial_parameters
    
    ############################################################################################################################################################################
    ############### choosing either the Baum-Welch algorithm or the method of directly maximizing the likelihood ###############################################################
    ############################################################################################################################################################################
    if (training_method == "EM") 
    {
      temp_trained_HMM_regarding_distribution_theta <- Baum_Welch_algorithm(x = x, m = temp_initial_parameters$m, gamma = temp_initial_parameters$gamma, delta = temp_initial_parameters$delta, distribution_class = distribution_class, distribution_theta = temp_initial_parameters$distribution_theta, discr_logL = discr_logL, Mstep_numerical = Mstep_numerical, discr_logL_eps = discr_logL_eps,  BW_max_iter = BW_max_iter, BW_limit_accuracy = BW_limit_accuracy, DNM_max_iter = DNM_max_iter, DNM_limit_accuracy = DNM_limit_accuracy,  DNM_print = DNM_print, BW_print = BW_print)
    }
    if (training_method == "numerical") 
    {
      temp_trained_HMM_regarding_distribution_theta = direct_numerical_maximization(x = x, m = temp_initial_parameters$m, distribution_class = distribution_class, gamma = temp_initial_parameters$gamma, distribution_theta = temp_initial_parameters$distribution_theta, DNM_max_iter = DNM_max_iter, DNM_limit_accuracy = DNM_limit_accuracy, DNM_print = DNM_print)
    }		
    list_of_all_trained_HMMs[[i]] <- temp_trained_HMM_regarding_distribution_theta
    
    list_of_all_logLs_for_each_HMM_with_m_states <- c(list_of_all_logLs_for_each_HMM_with_m_states, temp_trained_HMM_regarding_distribution_theta$logL)
    list_of_all_AICs_for_each_HMM_with_m_states <- c(list_of_all_AICs_for_each_HMM_with_m_states, temp_trained_HMM_regarding_distribution_theta$AIC)
    list_of_all_BICs_for_each_HMM_with_m_states <- c(list_of_all_BICs_for_each_HMM_with_m_states, temp_trained_HMM_regarding_distribution_theta$BIC)
  }	
  
  ############################################################################################################################################################################
  ############### Selection of the HMM with the most plausible m #############################################################################################################
  ############################################################################################################################################################################

  AIC_selection_of_m <- which.min(list_of_all_AICs_for_each_HMM_with_m_states)
  BIC_selection_of_m <- which.min(list_of_all_BICs_for_each_HMM_with_m_states)
  
  model_selection_over_AIC <- TRUE
  selected_m <- AIC_selection_of_m
  if(AIC_selection_of_m > BIC_selection_of_m)
  {
    selected_m <- BIC_selection_of_m
    model_selection_over_AIC <- FALSE
  }
  
  trained_HMM_with_selected_m <- list_of_all_trained_HMMs[[selected_m]]
  
  
  trained_HMM_with_selected_m <- list(x = trained_HMM_with_selected_m$x, 
                                      m = trained_HMM_with_selected_m$m, 
                                      logL = trained_HMM_with_selected_m$logL, 
                                      AIC = trained_HMM_with_selected_m$AIC, 
                                      BIC = trained_HMM_with_selected_m$BIC, 
                                      delta = trained_HMM_with_selected_m$delta, 
                                      gamma = trained_HMM_with_selected_m$gamma, 
                                      distribution_class = trained_HMM_with_selected_m$distribution_class, 
                                      distribution_theta= trained_HMM_with_selected_m$distribution_theta, 
                                      estimated_mean_values= trained_HMM_with_selected_m$estimated_mean_values)
                                      
############################################################################################################################################################################################################################### Return results ################################################################################################################
################################################################################################################################################################################
  
print(trained_HMM_with_selected_m)
  
return(list(trained_HMM_with_selected_m = trained_HMM_with_selected_m, 
            list_of_all_initial_parameters = list_of_all_initial_parameters, 
            list_of_all_trained_HMMs = list_of_all_trained_HMMs, 
            list_of_all_logLs_for_each_HMM_with_m_states = list_of_all_logLs_for_each_HMM_with_m_states, 
            list_of_all_AICs_for_each_HMM_with_m_states = list_of_all_AICs_for_each_HMM_with_m_states, 
            list_of_all_BICs_for_each_HMM_with_m_states = list_of_all_BICs_for_each_HMM_with_m_states, 
            model_selection_over_AIC=model_selection_over_AIC))
}
