#' The Three Biggest Drawdowns in the portfolio
#'
#' Show the top 3 drawdowns including start and end dates, as well as
#' decrease in returns during the drawdown period. All the information will be
#' returned in a table, with all the numbers properly formatted.
#' 
#' If the data set is not big enough that there are fewer drawdowns than 
#' required by the user, the function will throw NA's into the table so that the
#' table will still contain as many rows as the user demands.
#'
#' @param x A data frame with date and return columns
#' 
#' @return A data frame that contains the starting date, end date and values of
#' the three biggest drawdowns.

drawdown <- function(x) {
  
  ## Create a new column. Each value contains the biggest drawdowns that end 
  ## exactly at the corresponding date of that row
  
  res <- list()
  
  x[, "pnl"] <- cumsum(x[["pnl"]]) 
  
  diff <- cummax(x[["pnl"]]) - x[["pnl"]]
  ## Start the output table of the function with labels
  
  res[[1]] <- data.frame(start = "Start Date", end = "End Date", pnl = "P&L ($)")
  
  for (i in 1:3){
    
    ## Fill the output table with NA's if there is no more drawdown available
    
    if (length(diff) != 0){
      
      ## Find the end date of the biggest drawdown
      
      dd.end.idx   <- max(which.max(diff), length(diff) - which.max(rev(diff)) + 1)
    }
    if (nrow(x) > 1 & length(which(diff == 0 & seq(diff) < dd.end.idx)) != 0) {
      
      ## Look up the starting date of the biggest drawdown and calculate the 
      ## P&L loss of that drawdown
      
      dd.start.idx <- max(which(diff == 0 & seq(diff) < dd.end.idx))
      dd.pnl       <- x[["pnl"]][dd.end.idx] - x[["pnl"]][dd.start.idx]
      
      ## Find the dates of the recovery period from the biggest drawdown
      
      dd.recover.idx <- which(diff == 0 & seq(diff) > dd.end.idx)
      
      ## Find the start date of the recovery period
      
      if (length(dd.recover.idx) == 0) {
        dd.recover.idx <- length(diff)
      } else {
        dd.recover.idx <- min(dd.recover.idx)       
      }
      
      ## Store the results of the drawdown in a layer
      
      res[[i+1]] <- data.frame(start = as.character(next_trading_day(x[["date"]][dd.start.idx])),
                               end   = as.character(x[["date"]][dd.end.idx]),
                               pnl   = as.character(dd.pnl))
      
      ## Remove the data of the biggest drawdown and the recovery period
      
      x <- x[-1 * c(dd.start.idx:dd.recover.idx), ]
      diff <- diff[-1 * c(dd.start.idx:dd.recover.idx)]
      
    }
  }
  
  ## Make a table
  
  do.call("rbind", res)
  
}