#' Plot intraindividual variability of reaction time
#' 
#' Plot intraindividual variability in reaction time, faceted by the four essential blocks. 
#' 
#' @param my_data The raw dataframe to be used
#' @param session_id A string of the variable name identifying each unique participant.
#' @param data_type A string of "raw" for no cleaning, or "clean" for cleaned data (no error trials, RT < 10,000ms, and RT > 180ms)
#' @param block_name A string of the variable name for the blocks
#' @param trial_blocks A vector of the four essential blocks in the seven-block IAT (i.e., B3, B4, B6, and B7).
#' @param trial_number A string of the variable identifying the trial number.
#' @param trial_latency A string of the variable name for the latency of each trial.
#' @import ggplot2
#' @export

plotIIV <- function(my_data, data_type, block_name, trial_blocks, session_id, trial_number, trial_latency){
  
  if(length(unique(my_data[, session_id])) > 100){
    
    sample_ids <- sample(unique(my_data[, session_id]), 100)
    my_data <- my_data[my_data[, session_id] %in% sample_ids, ]
    warning("Your total sample size is > 100. A random subsample was taken for plotting.")
    
  }
  
  if(data_type == "clean"){
    
    data_new <- my_data[my_data[, block_name] %in% trial_blocks, ]
    
    data_new <- data_new[data_new[, trial_latency] > 180 & data_new[, trial_latency] < 10000, ]
    
  } else if(data_type == "raw"){
    
    data_new <- my_data[my_data[, block_name] %in% trial_blocks, ]
    
  } else stop("Please enter 'raw' or 'clean' for data_type.")
  
  data_new[, block_name] <- as.factor(data_new[, block_name])
  data_new[, trial_number] <- as.factor(data_new[, trial_number])
  
  ggplot(data_new, aes_string(x = trial_number, y = trial_latency)) + 
    geom_line(aes_string(group = session_id), color = "#999999") + 
    stat_smooth(method = "loess", aes(group = 1), size = 1.5, color = "black") +
    facet_grid(paste(". ~ ", block_name, sep = "")) + 
    theme(axis.text.x = element_blank()) + 
    theme_bw()

}

