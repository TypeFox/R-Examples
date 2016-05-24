rm(list = ls())

library(quantmod)
library(xts)

#' Function to simulate data for RMT
#' 
#' @param N number of observations
#' @param p number of variables
#' @param seed random seed
#' 
#' @author Rohit Arora
#' 
sim.RMTdata <- function(N = 500, p = 100, numOfSpikes = 15, seed = 9999) {
  
  if(p < numOfSpikes) 
    stop("Number of spikes must be less than number of variables")
  
  set.seed(seed)
  
  spikes <- seq(20, 48, length.out = numOfSpikes)
  E <- diag(spikes)

  M <- matrix(rnorm(N*p), nrow = N, ncol = p)
  C <- cov(M); vec <- eigen(C)$vectors[, 1:numOfSpikes]
  C <- vec %*% E %*% t(vec) + diag(1, p)
  
  data <- rmvnorm(N, sigma = C)
  
  to <- Sys.Date(); from <- Sys.Date() - N + 1
  dates <- seq(from = from, to = to, by = "days")
  xts(data, order.by = dates)
}



#' Function to download closing symbol data. This function must be used to download
#' and store data in the data folder when symbol list is updated or dates change. 
#' After executing the function store the symdata object in the data folder.
#' 
#' @details
#' This method takes in a list of symbols and downloads them for specific dates
#' @param symbols list of symbols in string format
#' @param from start date in yyyy-mm-dd format
#' @param to   end   date in yyyy-mm-dd format if missing then defaults to current date
#' @param returns boolean to indicate if prices must be converted to log returns
#' @param frequency string, "monthly"/"daily" to indicate the type of time series.
#'        defaults to monthly
#' @author Rohit Arora
#' 
get.closing.symdata <- function(symbols, from, to=Sys.Date(), returns=TRUE,
                                        frequency = c("monthly", "daily"))
  {
  
  if (length(symbols) == 0) stop("No symbols to download")
  
  d <- try( as.Date( from, format= "%Y-%m-%d" ) )
  if( class( d ) == "try-error" || is.na( d ) ) stop( "Invalid start date" )
  
  d <- try( as.Date( to, format= "%Y-%m-%d" ) )
  if( class( d ) == "try-error" || is.na( d ) ) stop( "Invalid end date" )  
  

  getSymbols(symbols, from = from, to = to)
  
  frequency = frequency[1]
  
  for(symbol in symbols) {
    
    if(symbol == "LOW") 
      colnames(LOW) <- gsub("^LOW\\.","",colnames(LOW))
    
    val <- get(symbol)
    
    temp0  <- temp <- adjustOHLC(val, symbol.name=symbol)
    if(frequency == "monthly") {
      temp    <- to.monthly(temp0, indexAt='endof', drop.time=FALSE)
      colnames(temp) <- colnames(temp0)      
    }
    assign(x=symbol, value=temp)
  }
  
  symbolData <- do.call(merge, lapply(symbols, function(x) Cl(get(x))))
  colnames(symbolData) <- symbols
  if (returns) symbolData <- diff(log(symbolData))[-1,]
  symbolData
}

menu <- function() {
  print("1 . Download symbols for factor data")  
  print("2 . Download symbols for filtering using RMT")  
  print("3 . Download symbols for Independent Switching DCC")  
  print("4 . Simulate data for Spiked Covariance Model")  
  print("5 . Exit")  
}

while(TRUE) {
  menu()
  opt <- readline("Select any option : ")
  
  switch(opt, 
         "1" = {
           symbols <-  c('BABA', 'TWTR', 'LNKD', 'YHOO', 'GE', 'LAZ', 'V')
           
           missingdata <- get.closing.symdata(symbols, from = "2007-04-01")
           l <- ls(); rm(list = c(l[l != "missingdata"],"l")) 
           break
         }, 
         "2" = {
           symbols <- c('AAPL','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE',
                        'GS','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM',
                        'MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','V',
                        'VZ','WMT','XOM')
           
           dow30data <- get.closing.symdata(symbols, from = "2014-04-01",
                                               frequency = "daily")
           dow30data <- na.omit(dow30data)
           l <- ls(); rm(list = c(l[l != "dow30data"],"l"))
           break
         },
         "3" = {
           symbols <- c('XLE', 'XLY', 'XLP','XLF','XLV','XLI','XLB','XLK', 'XLU')
           
           etfdata <- get.closing.symdata(symbols, from = "2008-01-01",
                                                to = "2010-12-31",
                                               frequency = "daily")
           etfdata <- na.omit(etfdata)
           l <- ls(); rm(list = c(l[l != "etfdata"],"l"))
           break
         },
         "4" = {
           rmtdata <- sim.RMTdata()
           l <- ls(); rm(list = c(l[l != "rmtdata"],"l"))
           break
         },
         "5" = {
            print('Exiting...')
           break
         },
         {
            print("Invalid Option. Please re-select")  
           cat("\n")
         } 
  )
}


