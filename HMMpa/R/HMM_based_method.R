HMM_based_method <-
function(x, cut_points, distribution_class, 
         min_m = 2, max_m = 6, n = 100,
         max_scaled_x = NA, names_activity_ranges = NA,  
         discr_logL = FALSE, discr_logL_eps = 0.5, 
         dynamical_selection = TRUE, training_method = "EM", 
         Mstep_numerical = FALSE, 
         BW_max_iter = 50, BW_limit_accuracy = 0.001, BW_print = TRUE,
         DNM_max_iter = 50, DNM_limit_accuracy = 0.001, DNM_print = 2, decoding_method = 'global',
         bout_lengths = NULL, plotting = 0)
{
	
################################################################################################################################################################################################################################# Check on arguments ##########################################################################################################
################################################################################################################################################################################  	
 if(is.null(bout_lengths))
 {
 	stop("Set variable 'bout_lengths' to use this function. See help-manual for further information. For example: bout_lengths=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,12,13,20,21,40,41,60,61,80,81,120,121,240,241,480,481,1440,1,1440)")
 }		   

################################################################################################################################################################################################################################# Needed variables and functions ##############################################################################################
################################################################################################################################################################################  	
  ############################################### function to scale counts in order to pretend from overflow/underflow-errors using training algorithms       #################

  scaling_observations <- function(x, max_scaled_x)
  {
    scaling_observations_factor <- max_scaled_x / max(x)
    
    scaled_x <- scaling_observations_factor * x
    
    return(list(original_x = x, scaling_observations_factor = scaling_observations_factor, scaled_x = scaled_x))
  }
  
  ############################################### scaling of the observations (counts)       ##################################################################################
  
  original_x <- x
  
  if(!is.na(max_scaled_x))
  {
    x <- scaling_observations(x = x, max_scaled_x = max_scaled_x)
    
    data_scale_factor <- x$scaling_observations_factor
    
    x <- x$scaled_x
  }
  if(distribution_class == "pois" | distribution_class == "genpois" | distribution_class == "bivariate_pois" | distribution_class == "geom")
  {
    x <- round(x)	
  }
  
################################################################################################################################################################################################################################# Training of (the most plausible) HMM for the given time-series of counts         ############################################
################################################################################################################################################################################
  
  trained_HMM_with_selected_m <- HMM_training(x = x, min_m = min_m, max_m = max_m, distribution_class = distribution_class, discr_logL = discr_logL, discr_logL_eps = discr_logL_eps, training_method = training_method, Mstep_numerical = Mstep_numerical, n = n, dynamical_selection = dynamical_selection,  BW_max_iter = BW_max_iter, BW_limit_accuracy = BW_limit_accuracy, BW_print = BW_print, DNM_max_iter = DNM_max_iter, DNM_limit_accuracy = DNM_limit_accuracy, DNM_print = DNM_print)$trained_HMM_with_selected_m

################################################################################################################################################################################################################################# Decoding ot the trained HMM for the given time-series of counts       #######################################################
################################################################################################################################################################################  
 

  	decoding <- HMM_decoding(x = x, m = trained_HMM_with_selected_m$m, delta = trained_HMM_with_selected_m$delta, gamma = trained_HMM_with_selected_m$gamma, distribution_class = trained_HMM_with_selected_m$distribution_class, distribution_theta = trained_HMM_with_selected_m$distribution_theta, decoding_method = decoding_method, discr_logL = discr_logL, discr_logL_eps = discr_logL_eps)
  
################################################################################################################################################################################################################################# back-scaling of the extracted hidden PA-levels underlying the time-series of counts       ####################################
################################################################################################################################################################################
    
  if(!is.na(max_scaled_x))
  {			
    decoding$decoding_distr_means <- (1 / data_scale_factor) * decoding$decoding_distr_means
  }else{
    decoding$decoding_distr_means <- decoding$decoding_distr_means
  }

################################################################################################################################################################################################################################# Applying the traditional cut-off point method on the extracted hidden PA-levels underlying the time-series of counts  #######
################################################################################################################################################################################
  
  extendend_cut_off_point_method <- cut_off_point_method(x = original_x, hidden_PA_levels = decoding$decoding_distr_means , cut_points = cut_points, names_activity_ranges = names_activity_ranges, bout_lengths = bout_lengths, plotting = plotting)
  
############################################################################################################################################################################################################################### Return results ################################################################################################################
################################################################################################################################################################################  

return(list(trained_HMM_with_selected_m = trained_HMM_with_selected_m,
            decoding = decoding,
            extendend_cut_off_point_method = extendend_cut_off_point_method))
  
}
