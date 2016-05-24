#' Find the best-performing and the worst-performing months
#' 
#' This function takes in a data set and returns a list with the best and the 
#' worst month and their respective pnl's.
#' 
#' @param x A data frame that contains data for individual instruments.
#'   
#' @return A list with the best month and the worst month, as well as their
#' respective p&l's.

best_worst_month <- function(x){
  
  ## For R CMD Check to not issue notes of visible binding stuff
  
  pnl    <- NULL
  
  ## Initiate the list to generate output
  
  output <- list()
  
  ## Select columns of the date and pnl in the data frame to save time for
  ## manipulating the data set. Then group the data set by year-month
  
  x <- x %>% 
    select(date, pnl) %>% 
    mutate(date = format(date, "%Y-%m")) %>%
    group_by(date) %>%
    summarise(pnl = sum(pnl, na.rm = TRUE)) %>%
    ungroup()
  
  ## Find the year-month with the largest p&l. Then find the pnl of the best month
  
  output$best.month <- x$date[which.max(x$pnl)]
  output$best.pnl <- max(x$pnl)
  
  ## Find the worst month and its p&l
  
  output$worst.month <- x$date[which.min(x$pnl)]
  output$worst.pnl <- min(x$pnl)
  
  return(output)
}
