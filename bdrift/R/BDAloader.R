#' @title Beta Drift Anaylsis Data Loader
#'
#' @description
#' \code{BDA.loader} prepares a data frame to be used 
#' by the \code{BDA} function.
#'
#' @details
#' \code{BDA.loader} pulls stock price data from Yahoo Finance, 
#' calculates returns on these prices, downloads factor data from 
#' Kenneth French's library (via Quandl.com) and bundles all 
#' data in a xts-matrix that can be passed on to the \code{BDA} 
#' function.
#' 
#' 
#' @param symbol stock ticker on Yahoo Finance, 
#' enter as character.
#' @param frequency the frequency used to calculate returns 
#' (\code{"daily"}, \code{"monthly"}, or \code{"yearly"})
#' @param xbench ticker symbol of an external benchmark, 
#' \code{NA} by default.
#' @param type type of returns to be calculated (\code{"log"} 
#' or \code{"arithmetic"}). By default, log returns are used.
#' @param ... aditional commands passed to the \code{getSymbols} 
#' function.
#' @export
#' @return a xts-matrix containing the returns of the security, 
#' Kenneth French's asset pricing factors and the external 
#' benchmark (optional).
#' 
#' @author Markus Peter Auer <mp.auer@@meanerreversion.com>
#' @examples 
#' testframe <- BDA.loader(symbol = "XOM")

BDA.loader <- function(symbol, 
                        frequency = "monthly",
                        xbench = NA,
                        type = "log",
                        ...){
  if (frequency == "daily") {
    FF <- Quandl::Quandl("KFRENCH/FACTORS5_D")
    MOM <- Quandl::Quandl("KFRENCH/MOMENTUM_D")
    FF$Date <- zoo::as.Date(as.character(FF$Date), format='%Y-%m-%d')
    MOM$Date <- zoo::as.Date(as.character(MOM$Date), format='%Y-%m-%d')
    FF <- xts::as.xts(FF[,2:7], order.by = FF$Date)/100
    MOM <- xts::as.xts(MOM[,2], order.by = MOM$Date)/100
    factors <- xts::merge.xts(FF, MOM, join = "inner")
  } else if (frequency == "monthly") {
    FF <- Quandl::Quandl("KFRENCH/FACTORS5_M")
    MOM <- Quandl::Quandl("KFRENCH/MOMENTUM_M")
    FF$Date <- zoo::as.Date(as.character(FF$Date), format='%Y-%m-%d')
    MOM$Date <- zoo::as.Date(as.character(MOM$Date), format='%Y-%m-%d')
    FF <- xts::as.xts(FF[,2:7], order.by = FF$Date)/100
    MOM <- xts::as.xts(MOM[,2], order.by = MOM$Date)/100
    factors <- xts::merge.xts(FF, MOM, join = "inner")
    factors <- xts::to.monthly(factors, OHLC = FALSE)
  } else if (frequency == "yearly") {
    FF <- Quandl::Quandl("KFRENCH/FACTORS5_A")
    MOM <- Quandl::Quandl("KFRENCH/MOMENTUM_A")
    FF$Date <- zoo::as.Date(as.character(FF$Date), format='%Y-%m-%d')
    MOM$Date <- zoo::as.Date(as.character(MOM$Date), format='%Y-%m-%d')
    FF <- xts::as.xts(FF[,2:7], order.by = FF$Date)/100
    MOM <- xts::as.xts(MOM[,2], order.by = MOM$Date)/100
    factors <- xts::merge.xts(FF, MOM, join = "inner")
    factors <- xts::to.yearly(factors, OHLC = FALSE)
  }
  asset <- quantmod::Ad(quantmod::getSymbols(Symbols = symbol, auto.assign = FALSE, ...))
  cna <- gsub(".Adjusted", "", colnames(asset))
  asset <- quantmod::periodReturn(x = asset, 
                        period = frequency, 
                        type = type) 
  colnames(asset) <- cna
  symbol <- cna
  
  if (is.na(xbench) == FALSE) { 
    bench <- quantmod::Ad(quantmod::getSymbols(Symbols = xbench, auto.assign = FALSE, ...))
    cnb <- gsub(".Adjusted", "", colnames(bench))
    bench <- quantmod::periodReturn(x = bench, 
                          period = frequency, 
                          type = type)
    colnames(bench) <- cnb
    xbench <- cnb
    asset <- xts::merge.xts(asset, bench, join = "inner")
  }
  if (frequency == "monthly") {
    asset <- xts::to.monthly(asset, OHLC = FALSE)
  } else if (frequency == "yearly") {
    asset <- xts::to.yearly(asset, OHLC = FALSE)
  }
  DF <- xts::merge.xts(factors, asset, join = "inner")
  DF[,symbol] <- DF[,symbol]-DF$RF
  if (is.na(xbench) == FALSE) {
    DF[,xbench] <- DF[,xbench]-DF$RF
  }
  return(DF)
}