#' Summarize data for the overall portfolio
#' 
#' This function takes in a data set and returns the summary for each day's
#' net market value, profit and loss, contract number and return rate.
#' 
#' @param x A data frame that contains data for individual commodities.
#' 
#' @param sector.selected The sector that the function is doing summarization
#' for. When "sector.selected = NULL", the function performs summarization across
#' the whole data set.
#' 
#' @param capital.num The constant number of allocated capital to calculate return, 
#' number. Default value is NULL. When "capital.num = NULL", the function uses
#' gross market value to calculate return rates for each row.
#' 
#' @return A data frame that summarizes the market values and returns across all
#' commodities.

sum_data <- function(x, sector.selected = NULL, capital.num = NULL){
  
  ## initialization so that CMD Check won't release notes
  
  sector <- nmv <- pnl <- contract <- gmv <- NULL
  
  ## Select a specific sector if this function is asked to summarize across a
  ## a specific sector
  
  if(!is.null(sector.selected)) {
    x <- filter(x, sector == sector.selected)
  }
  
  if(dim(x)[1] == 0){
    stop("There does not exist any data that satisfies your selection. 
         Please select a new combination of Strategy, Portfolio and Instrument.")
  }
  
  ## Group the data set by date 
  
  x <- group_by(x, date)
  
  ## Calculate the total nmv, gmv, pnl, and contract for the data set
  
  x.sum <- summarise(x,
                     nmv = sum(nmv, na.rm = TRUE),
                     pnl = sum(pnl, na.rm = TRUE),
                     contract = sum(contract, na.rm=TRUE),
                     gmv = sum(gmv, na.rm = TRUE))
  x.sum <- tbl_df(x.sum) %>%
    arrange(date) %>%
    filter(gmv != 0)
  
  ## Calculate return rates with allcated capital if the user has that variable,
  ## or calculate with gross market value if capital.num is unavailable
  
  if(is.null(capital.num)){
    x.sum <- mutate(x.sum, ret = pnl / gmv, cumpnl = cumsum(pnl))
  } else {
    x.sum <- mutate(x.sum, ret = pnl / capital.num, cumpnl = cumsum(pnl))
  }
  
  ## Return the new data frame with new column names and values
  
  x.sum <- mutate(x.sum, name = " ")
  
  ## Generate a name entry so that the graph functions can take the name and
  ## draw a graph
  
  if(!is.null(sector.selected)){
    x.sum[[1,8]] <- sector.selected
  } else {
    x.sum[[1,8]] <- "the Whole Portfolio"
  }
  
  return(x.sum)
  
}
