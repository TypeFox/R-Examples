#Delay of up to 1 minute returns all markets

GeneralMarketDataAll <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api.php?method=marketdatav2 ")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' DOGECOIN to USD
#'
#' This function allows you to get general market data on Dogecoin to US Dollars
#' @keywords doge
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_DOGEUSD()
#' }


#DOGECOIN to USD - Realtime

GeneralMarketData_DOGEUSD <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/182")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' BITCOIN to USD
#'
#' This function allows you to get general market data on Bitcoin to US Dollars
#' @keywords bitcoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_BTCUSD()
#' }

#BITCOIN to USD - Realtime

GeneralMarketData_BTCUSD <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/2")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Feathercoin to USD
#'
#' This function allows you to get general market data on Feathercoin to US Dollars
#' @keywords feathercoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_FTCUSD()
#' }

#FEATHERCOIN to USD - Realtime

GeneralMarketData_FTCUSD <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/6")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Litecoin to USD
#'
#' This function allows you to get general market data on Litecoin to US Dollars
#' @keywords litecoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_LTCUSD()
#' }

#LITECOIN to USD - Realtime

GeneralMarketData_LTCUSD <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/1")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Darkcoin to USD
#'
#' This function allows you to get general market data on Darkcoin to US Dollars
#' @keywords Darkcoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_DRKUSD()
#' }

#DARKCOIN to USD - Realtime

GeneralMarketData_DRKUSD <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/213")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Ripple to USD
#'
#' This function allows you to get general market data on Ripple to US Dollars
#' @keywords Ripple
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_XRPUSD()
#' }

#Ripple to US Dollars- Realtime

GeneralMarketData_XRPUSD <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/442")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
  
}

#' ReddCoin to USD
#'
#' This function allows you to get general market data on ReddCoin to US Dollars
#' @keywords ReddCoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_RDDUSD()
#' }

#ReddCoin to US Dollars- Realtime

GeneralMarketData_RDDUSD <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/262")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
  
}

#' Peercoin to USD
#'
#' This function allows you to get general market data on Peercoin to US Dollars
#' @keywords Peercoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_PCCUSD()
#' }

#PeerCoin to US Dollars- Realtime

GeneralMarketData_PCCUSD <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/305")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
  
}

#' Dogecoin to Bitcoin
#'
#' This function allows you to get general market data on Dogecoin to BITCOIN
#' @keywords Dogecoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_DOGEBTC()
#' }

#DOGECOIN to BITCOIN - Realtime

GeneralMarketData_DOGEBTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/132")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Darkcoin to Bitcoin
#'
#' This function allows you to get general market data on Darkcoin to BITCOIN
#' @keywords Darkcoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_DRKBTC()
#' }

#DARKCOIN to BITCOIN - Realtime

GeneralMarketData_DRKBTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/2")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Feathercoin to Bitcoin
#'
#' This function allows you to get general market data on Feathercoin to BITCOIN
#' @keywords Feathercoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_FTCBTC()
#' }

#FEATHERCOIN to BITCOIN - Realtime

GeneralMarketData_FTCBTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/5")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Litecoin to Bitcoin
#'
#' This function allows you to get general market data on Litecoin to BITCOIN
#' @keywords Litecoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_LTCBTC()
#' }

#LITECOIN to BITCOIN - Realtime

GeneralMarketData_LTCBTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/3")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' NXTcoin to Bitcoin
#'
#' This function allows you to get general market data on NXTcoin to BITCOIN
#' @keywords NXTcoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_NXTBTC()
#' }

#NXTCOIN to BITCOIN - Realtime

GeneralMarketData_NXTBTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/159")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' 42coin to Bitcoin
#'
#' This function allows you to get general market data on 42coin to BITCOIN
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_42CBTC()
#' }

#42COIN to BITCOIN - Realtime

GeneralMarketData_42CBTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/141")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Unobtanium to Bitcoin 
#'
#' This function allows you to get general market data on Unobtaniumto US Dollars
#' @keywords Unobtanium
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_UNOBTC()
#' }


#Unobtanium to Bitcoin - Realtime

GeneralMarketData_UNOBTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/133")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' Dogecoin to Litecoin
#'
#' This function allows you to get general market data on Dogecoin to LITECOIN
#' @keywords Dogecoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_DOGELTC()
#' }

#DOGECOIN to LITECOIN - Realtime

GeneralMarketData_DOGELTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/135")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' AndroidsTokensV2  to Litecoin
#'
#' This function allows you to get general market data on AndroidsTokensV2 to LITECOIN
#' @keywords AndroidsTokensV2 
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_ADTLTC()
#' }

#AndroidsTokensV2 to LITECOIN - Realtime

GeneralMarketData_ADTLTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/94")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' AnonCoin  to Litecoin
#'
#' This function allows you to get general market data on AnonCoin to LITECOIN
#' @keywords AnonCoin 
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_ANCLTC()
#' }

#AnonCoin to LITECOIN - Realtime

GeneralMarketData_ANCLTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/121")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' AsicCoin to Litecoin
#'
#' This function allows you to get general market data on AsicCoin to LITECOIN
#' @keywords AsicCoin  
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_ASCLTC()
#' }

#AsicCoin  to LITECOIN - Realtime

GeneralMarketData_ASCLTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/111")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' AuroraCoin to Litecoin
#'
#' This function allows you to get general market data on AuroraCoin to LITECOIN
#' @keywords AuroraCoin   
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_AURLTC()
#' }

#AuroraCoin  to LITECOIN - Realtime

GeneralMarketData_AURLTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/161")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' BatCoin to Litecoin
#'
#' This function allows you to get general market data on BatCoin to LITECOIN
#' @keywords BatCoin  
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_BATLTC()
#' }

#BatCoin to LITECOIN - Realtime

GeneralMarketData_BATLTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/186")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' BlackCoin to Litecoin
#'
#' This function allows you to get general market data on BlackCoin  to LITECOIN
#' @keywords BlackCoin
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_BATLTC()
#' }

#BlackCoin  to LITECOIN - Realtime

GeneralMarketData_BCLTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/191")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}


#' CAIx to Litecoin
#'
#' This function allows you to get general market data on CAIx to LITECOIN
#' @keywords CAIx
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_CAIxLTC()
#' }

#CAIx  to LITECOIN - Realtime

GeneralMarketData_CAIxLTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/222")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}

#' CryptogenicBullion to Litecoin
#'
#' This function allows you to get general market data on CryptogenicBullion to LITECOIN
#' @keywords CryptogenicBullion
#' @export
#' @examples
#' \dontrun{
#' GeneralMarketData_CGBLTC()
#' }

#CryptogenicBullion  to LITECOIN - Realtime

GeneralMarketData_CGBLTC <- function () {
  internetcheck <- url.exists("https://api.cryptsy.com", timeout = 30)
  if( internetcheck != TRUE)
    stop('Cryptsy or your internet connection is down')
  data <- getURL("https://api.cryptsy.com/api/v2/markets/123")
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}
