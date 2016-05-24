#' Plot individual variability in the IAT
#' 
#' Plot mean participant reaction time with 95\% confidence intervals to see how reaction time varies by participant. The data is automatically cleaned
#' (i.e., no error trials, trials with RT > 10000 or < 180 are deleted) to avoid over/underinflation of mean estimates and only includes trials from
#' essential blocks.
#' 
#' @param my_data The raw dataframe to be used
#' @param session_id A string of the variable name identifying each unique participant.
#' @param block_name A string of the variable name for the blocks
#' @param trial_blocks A vector of the four essential blocks in the seven-block IAT (i.e., B3, B4, B6, and B7).
#' @param trial_latency A string of the variable name for the latency of each trial.
#' @param trial_error A string of the variable name identifying whether a trial was an error or not (1 = error)
#' @import ggplot2 dplyr
#' @importFrom stats reshape sd
#' @export

plotIndVar <- function(my_data, block_name, trial_blocks, session_id, trial_latency, trial_error){
  
  if(length(unique(my_data[, session_id])) > 100){
    
    sample_ids <- sample(unique(my_data[, session_id]), 100)
    
    my_data <- my_data %>%
      filter_(interp(~ x %in% sample_ids, x = as.name(session_id)))
    
    warning("Your total sample size is > 100. A random subsample was taken for plotting.")
    
  }
  
  by_part <- my_data %>%
    group_by_(session_id) %>%
    filter_(interp(~ x > 180 & x < 10000, x = as.name(trial_latency)),
            interp(~ y == 0, y = as.name(trial_error)),
            interp(~ z %in% trial_blocks, z = as.name(block_name))) %>%
    summarise_(mean_rt = interp(~ mean(x), x = as.name(trial_latency)),
               rt_se = interp(~ sd(x)/length(x), x = as.name(trial_latency))) %>%
    mutate_(session_id = interp(~ factor(a, levels = a[order(mean_rt)]),
                                a = as.name(session_id)))
  
  ggplot(by_part, aes_string(x = "session_id", y = "mean_rt")) + 
    geom_errorbar(aes(ymin =  by_part$mean_rt - 2 * by_part$rt_se,
                      ymax = by_part$mean_rt + 2 * by_part$rt_se), width = 0) +
    geom_point() +
    labs(x = "Session ID\n", y = "\nMean Reaction Time") + 
    coord_flip() + 
    theme_bw()
  
}
