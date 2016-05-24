#' Plot proportion of errors per item in the IAT
#' 
#' Plot proportion of errors in the IAT to see if any items yield high error rates. The data is automatically cleaned
#' (i.e., trials with RT > 10000 or < 180 are deleted) to avoid over/underinflation of mean error estimates
#' 
#' @param my_data The raw dataframe to be used
#' @param item_name A string of the variable identifying the items
#' @param trial_latency A string of the variable name for the latency of each trial.
#' @param trial_error A string of the variable name identifying whether a trial was an error or not (1 = error)
#' @import ggplot2 dplyr lazyeval
#' @export

plotItemErr <- function(my_data, item_name, trial_latency, trial_error){
  
  by_item_error <- my_data %>%
    group_by_(item_name) %>%
    filter_(interp(~ x > 180 & x < 10000, x = as.name(trial_latency))) %>%
    summarise_(prop_errors = interp(~ sum(y == 1)/length(y), y = as.name(trial_error))) %>%
    mutate_(item_name = interp(~ factor(z, levels = z[order(prop_errors)]),
                               z = as.name(item_name)))
  
  ggplot(by_item_error, aes_string(x = "item_name", y = "prop_errors")) + 
    geom_bar(fill = "dark gray", stat = "identity") +
    geom_hline(aes(yintercept = mean(by_item_error$prop_errors)), size = 1.5) +
    labs(x = "Item\n", y = "\nProportion Errors") + 
    coord_flip() + 
    theme_bw()
    
}