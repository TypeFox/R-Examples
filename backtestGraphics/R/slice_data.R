#' Slice Data Set
#' 
#' This function is part of the Shiny reactive functions. The function looks at
#' what strategy, portfolio and instrument the user selects, and then slices 
#' the data set into different pieces or summarizes the data set into a summary.
#' 
#' @param x is an input data frame that is to be sliced according to the user's
#' selection of strategy, portfolio and instrument.
#' @param input is a Shiny object that contains all the shiny inputs as a list.
#' @param capital.num is the amount of allocated capital for the portfolio, 
#' numeric type. The function will use capital.num to calculate return rates
#' if available, or use everyday's GMV to calculate return rates instead.
#' 
#' @return A small data set that has been retrieved from the previous data set.

slice_data <- function(x, input, capital.num){
  
  ## Set up vatiables to fool R CMD check
  
  name <- id <- portfolio <- gmv <- nmv <- contract <- x.name <- strategy <- NULL
  substrategy <- capital.num <- pnl <- sector <- NULL
  
  ## Generate a data frame of names that helps later selections
  
  x.name <- unique(select(x, name, id, sector, portfolio, strategy, substrategy))
  
  f <- list()
  
  ## Filter the correct strategy and portfolio, as specified by the user in 
  ## Shiny app. The conditions observe what the user selects in the three
  ## dropdown menus. For the strategy part, the function first determines if
  ## the user selects "strategy summary", a strategy or a substrategy, and then 
  ## either summary across all strategies for "strategy summary" or filter the
  ## required portion of strategy/substrategy. Then the function observes
  ## portfolios and instruments and summarize/filter the data frame according to
  ## these conditions.
  
  if(input$strategy == "Strategy Summary"){
    f$x.temp <- x %>% 
      group_by(name, id, date, sector, portfolio) %>%
      summarise(gmv = sum(gmv, na.rm = TRUE),
                nmv = sum(nmv, na.rm = TRUE),
                pnl = sum(pnl, na.rm = TRUE),
                contract = sum(abs(contract), na.rm = TRUE)) %>% 
      ungroup()
  } 
  else if (input$strategy %in% x.name$strategy){
    f$x.temp <- x %>% 
      filter(strategy == input$strategy)
  } 
  else {
    f$x.temp <- x %>% 
      filter(substrategy == input$strategy)
  }
  
  if(input$portfolio == "Portfolio Summary"){
    f$x.temp <- f$x.temp %>% 
      group_by(name, id, date, sector) %>%
      summarise(gmv = sum(gmv, na.rm = TRUE),
                nmv = sum(nmv, na.rm = TRUE),
                pnl = sum(pnl, na.rm = TRUE),
                contract = sum(abs(contract), na.rm = TRUE)) %>% 
      ungroup()
  } 
  else {
    f$x.temp <- f$x.temp %>% 
      filter(portfolio == input$portfolio)
  }
  
  if(input$instrument == "Instrument Summary"){
    f$x <- sum_data(f$x.temp,
                    capital.num = capital.num)
    
    f$instrument <- length(unique(f$x.temp[["id"]]))
  }
  
  if(input$instrument %in% unique(x.name$sector)){
    
    f$x <- sum_data(x = f$x.temp, 
                    sector.selected = input$instrument,
                    capital.num = capital.num)
    
    ## Use the original big data set to extract the 
    ## number of instruments for the sector because the sum_data function does
    ## not return instrument names after the summary
    
    f$sect <- filter(f$x.temp, sector == input$instrument)
    f$instrument <- length(unique(f$sect[["id"]]))
  }
  
  ## Do the same process for individual commodities. If the user has
  ## allocated capital, use that to calculate return rates. Otherwise
  ## use gross market value
  
  if(input$instrument %in% unique(x.name$id)){
    if(!is.null(capital.num)){
      f$x <- filter(f$x.temp, id == input$instrument) %>%
        filter(gmv != 0) %>% 
        mutate(cumpnl = cumsum(pnl), ret = pnl / capital.num) %>%
        as.data.frame()
    } else {
      f$x <- filter(f$x.temp, id == input$instrument) %>%
        filter(gmv != 0) %>% 
        mutate(cumpnl = cumsum(pnl), ret = pnl / gmv) %>%
        as.data.frame()
    }
    f$instrument <- 1
  }
  
  ## Check if the filtered data is empty. and return an error if it is
  
  validate(
    need(is.data.frame(f$x) | dim(f$x)[1] != 0, 
         "There does not exist any data that satisfies your selection. 
         Please select a new combination of Strategy, Portfolio and Instrument.")
  )
  
  return(f)
  
}