################################################################################
##
## $Id: calc.turnover.R 1300 2008-08-27 21:01:11Z zhao $
##
## Calculates the turnover for a long-short portfolio
## corresponding to a backtest.
##
################################################################################

## NOTE: This measure of turnover ignores price changes

## "x" is a character vector of unique identifiers

## "portfolio.factor" is a factor created from a vector of numeric
##  values as specified by an in.var

## "date.factor" is a factor created from a vector of dates as
##  specified by date.var

calc.turnover <- function(x, portfolio.factor, date.factor){
  
  stopifnot(
            is.factor(portfolio.factor),
            is.factor(date.factor),
            all.equal(length(x), length(portfolio.factor)),
            all.equal(length(x), length(date.factor))
            )
  
  data <- data.frame(x, portfolio.factor, date.factor)
  
  ## Split data by buckets
  
  data <- split(data, portfolio.factor)
  data <- data[sort(names(data))]
  
  data[["high"]]$longshort <- "long"
  data[["low"]]$longshort <- "short"
  
  ## Discard all rows except those in the "high" and "low" buckets
  
  data <- rbind(data[["high"]], data[["low"]])
  
  data$longshort <- factor(data$longshort, levels = c("long", "short"))
  
  ## Split by date
  
  data <- split(data, data$date.factor)
  data <- data[sort(names(data))]
  
  ## Array of turnovers for each date
  
  result <- array(dim = length(data), dimnames = list(names(data)))
  
  ## First date has no turnover
  
  last.portfolio <- list(subset(data[[1]], data[[1]]$longshort == "long"),
                         subset(data[[1]], data[[1]]$longshort == "short"))
  
  result[1] <- NA
  
  ## Loop through dates
  
  if(length(data) > 1){
    for(k in 2:length(data)){
      
      ## Split next portfolio by long/short
      
      new.portfolio <- list(subset(data[[k]], data[[k]]$longshort == "long"),
                            subset(data[[k]], data[[k]]$longshort == "short"))
      
      turnovers <- numeric()
      
      ## Loop both sides of portfolio
      
      for(i in 1:2){
        
        ## Which stocks are in this side of the portfolio?
        
        stocks <- union(as.character(new.portfolio[[i]]$x),
                        as.character(last.portfolio[[i]]$x))
        
        ## Create stock table for this side of the portfolio
        
        stock.table <- array(0, dim = c(length(stocks), 2),
                             dimnames = list(stocks, c("last", "new")))
        
        ## Calculate weights for last portfolio
        
        stock.table[as.character(last.portfolio[[i]]$x),"last"] <-
          100/nrow(last.portfolio[[i]])
        
        ## Calculate weights for new portfolio
        
        stock.table[as.character(new.portfolio[[i]]$x),"new"] <-
          100/nrow(new.portfolio[[i]])
        
        ## Store amount traded for each stock
        
        turnovers <- c(turnovers, abs(stock.table[ ,"last"] -
                                      stock.table[ ,"new"]))
      }
      
      ## Calculate turnover
      
      result[k] <- sum(turnovers)/400
      
      last.portfolio <- new.portfolio
      
    }
  }  
  invisible(result)
  
}
