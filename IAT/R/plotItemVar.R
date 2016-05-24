#' Plot IAT item variability
#' 
#' Plot mean item reaction time with 95\% confidence intervals to see how reaction time varies by item. The data is automatically cleaned
#' (i.e., no error trials, and trials with RT > 10000 or < 180 are deleted) to avoid over/underinflation of mean estimates and only includes
#'trials from essential blocks.
#' 
#' @param my_data The raw dataframe to be used
#' @param block_name A string of the variable name for the blocks
#' @param trial_blocks A vector of the four essential blocks in the seven-block IAT (i.e., B3, B4, B6, and B7).
#' @param item_name A string of the variable identifying the items
#' @param trial_latency A string of the variable name for the latency of each trial.
#' @param trial_error A string of the variable name identifying whether a trial was an error or not (1 = error)
#' @import ggplot2 dplyr
#' @importFrom stats reshape sd
#' @export

plotItemVar <- function(my_data, block_name, trial_blocks, item_name, trial_latency, trial_error){
  
  by_item_var <- my_data %>%
    group_by_(item_name) %>%
    filter_(interp(~ x > 180 & x < 10000, x = as.name(trial_latency)),
            interp(~ y == 0, y = as.name(trial_error)),
            interp(~ z %in% trial_blocks, z = as.name(block_name))) %>%
    summarise_(mean_rt = interp(~ mean(x), x = as.name(trial_latency)),
               rt_se = interp(~ sd(x)/length(x), x = as.name(trial_latency))) %>%
    mutate_(item_name = interp(~ factor(a, levels = a[order(mean_rt)]),
                               a = as.name(item_name)))
  
  ggplot(by_item_var, aes_string(x = "item_name", y = "mean_rt")) +
    geom_errorbar(aes(ymin = by_item_var$mean_rt - 2 * by_item_var$rt_se,
                      ymax = by_item_var$mean_rt + 2 * by_item_var$rt_se), width = 0) +
    geom_point() +
    labs(x = "\nItem", y = "Mean Reaction Time\n") +
    coord_flip() +
    theme_bw()
  
}
