#' TFX package: An \R interface to the TrueFX(tm) Market Data Web API 
#' 
#' \emph{TrueFX(tm) is a brand name that is owned by Integral Development Corp.
#' This software is in no way affiliated, endorsed, or approved by 
#' TrueFX(tm), Integral Development Corp, or any of its affiliates.}
#' The \pkg{TFX} package provides functions to access the free, real-time,
#' top-of-book, dealable, interbank foreign exchange rate quotes provided by 
#' TrueFX(tm).
#' 
#' The Web API can be accessed with an authenticated or an unauthenticated 
#' session.  The following is from TrueFX(tm) Market Data Web API Developer 
#' Guide: "Unauthenticated requests ... return a default response in the
#' form of a snapshot session that includes data for all currency pairs...
#' Unauthenticated requests do not support the more powerful query-based
#' functionality, including incremental responses containing only the changed
#' currency pairs."
#' 
#' A registered TrueFX(tm) account with a confirmed username is required to 
#' create an authenticated session (TrueFX(tm) accounts are free).  
#'
#' There are typically three steps to creating an authenticated TrueFX(tm)
#' session and requesting data:
#' 
#' \code{\link{ConnectTrueFX}} is used to request a \code{TFXsession} object 
#' which has a session ID needed to make a data request.
#'
#' \code{\link{QueryTrueFX}} will request market data from TrueFX(tm) using a 
#' \code{TFXsession} object that was created by \code{ConnectTrueFX}
#' 
#' \code{\link{ParseTrueFX}} will parse the results returned by 
#' \code{QueryTrueFX} into something that is easier to work with.
#'
#' There are also functions to \code{\link{Disconnect}} and 
#' \code{\link{Reconnect}} a \code{TFXsession}
#'
#' This package does not yet have explicit support streaming data, but it can be
#' accomplished a few different ways, one of which would be to use a 
#' \code{while} loop.
#' 
#' In addition to real time data, TrueFX(tm) also offers historical data since 
#' 2009.
#' 
#' From the TrueFX(tm) website,
#' "TrueFX is the first service that brings you real, dealable prices from real 
#' market participants from all the major market makers, with absolutely no 
#' intermediary. As a technology company, we can offer you historical 
#' tick-by-tick market data, at zero cost to you.
#' 
#' This data is top-of-the-book, tick-by-tick market data, with fractional pip 
#' spreads in millisecond detail. All timestamps are based on GMT."
#' 
#' The data can be downloaded from \url{http://www.truefx.com/?page=downloads}
#' 
#' There is a script in the /inst/parser directory of the FinancialInstrument
#' pacakage (\url{http://tinyurl.com/DownloadTrueFX}) that will download all 
#' available data, convert it to \code{xts} and save it to disk in binary files 
#' that are split by day such that \code{FinancialInstrument:::getSymbols.FI} 
#' can easily read them.
#' 
#' Some version of that script may make its way into a future release of this 
#' package.
#' 
#' @name TFX
#' @aliases TFX TFX-package
#' @docType package
#' @author Garrett See \email{gsee000@@gmail.com}
#' @references 
#' \url{http://www.truefx.com}
#'
#' \url{http://www.truefx.com/dev/data/TrueFX_MarketDataWebAPI_DeveloperGuide.pdf}
#' @keywords package programming IO
#' @examples
#' \dontrun{
#' ## Unauthenticated
#' QueryTrueFX()
#' QueryTrueFX(pretty=FALSE)
#' 
#' ## Must have a TrueFX(tm) account to run the following (Membership is free)
#' ## Replace JSTrader and Ou812 with your username and password, respectively
#' id <- ConnectTrueFX("EUR/USD,AUD/JPY", u='JSTrader', p='Ou812', f='html')
#' QueryTrueFX(id)
#' 
#' browseURL(paste0("http://webrates.truefx.com/rates/connect.html?id=", id))
#'
#' #view the Web API Developer Guide:
#' TrueFXRef()
#' }
NA
