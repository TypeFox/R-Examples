#' Find the three best-performing and worst-performing commodities
#' 
#' This function takes in a data set and returns the best three and worst three 
#' commodity names and their respective pnls. All these data will be presented
#' as a formatted table.
#'
#' @param x A data frame that contains data for individual commodities.
#' 
#' @return output A list of the best three and worst three performers with their 
#' commodity names, ids and respective pnls

best_worst_three <- function(x){
  
  ## Initiate so that CMD Check won't release notes
  
  name <- pnl <- desc <- NULL
  
  ## Initiate a list to contain all outputs, and the first rows of the two output
  ## data frames
  
  output <- best.frame <- worst.frame <- list()
  best.frame[[1]]  <- data.frame(Instruments = "Instruments", pnl = "P&L ($)")
  worst.frame[[1]] <- data.frame(Instruments = "Instruments", pnl = "P&L ($)")
  
  ## Select columns of the name and pnl in the data frame. Then group the data
  ## set by name, summarize pnl for individual instruments, and order instruments
  ## by pnl value
  
  x <- x %>% 
    select(pnl, name) %>%
    group_by(name) %>%
    summarise(pnl = sum(pnl, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(pnl))
  
  ## Take out the best/worst performers, then the second best/worst and the
  ## third. If there are less than three instruments available, the function
  ## will only look at all the available ones
  
  n <- nrow(x)
  for(i in 1:min(3, length(unique(x$name)))){
    best.frame[[i+1]]  <- data.frame(Instruments = x$name[i], 
                                     pnl = as.character(x$pnl[i]))
    worst.frame[[i+1]] <- data.frame(Instruments = x$name[n - i + 1], 
                                     pnl = as.character(x$pnl[n - i + 1]))
  }
  
  output$best.3  <- do.call("rbind", best.frame)
  output$worst.3 <- do.call("rbind", worst.frame)
  
  return(output)   
}