
#' @export
#' @rdname QueryTrueFX
ConnectTrueFX <- function(currency.pairs, username, password, 
                          qualifier='default', format, snapshot=FALSE) {
  if (missing(format)) format <- "default"
  if (!substr(format[[1L]], 1, 1) %in% c("d", "c", "h")) {
    warning("unrecognized format. Using default")
    format <- 'default'
  }
  if (missing(currency.pairs) || 
      (!missing(currency.pairs) && nchar(currency.pairs) < 1)) {
    ## If missing, use the 15 pairs for which TrueFX(tm) offers historical data
    currency.pairs <- c("EUR/USD", "USD/JPY", "GBP/USD", "EUR/GBP", "USD/CHF", 
                        "EUR/JPY", "EUR/CHF", "USD/CAD", "AUD/USD", "GBP/JPY", 
                        "AUD/JPY", "AUD/NZD", "CAD/JPY", "CHF/JPY", "NZD/USD")
    ## However, an unauthenticated request only returns the first 10 of those
  }
  stopifnot(is.character(currency.pairs))
  x <- unlist(strsplit(gsub(" ", "", currency.pairs), ","))
  cp <- paste(paste(substring(x, 1, 3), 
                    substring(x, nchar(x)-2, nchar(x)), sep="/"), 
              collapse=",")
  base.url <- "http://webrates.truefx.com/rates/connect.html"
  session <- new.env()
  
  if (missing(username) || missing(password)) {
    stop("missing username or password.")
  } else {
    URL <- paste0(base.url, 
                  "?u=", username, 
                  "&p=", password, 
                  "&q=", qualifier, 
                  "&c=", cp)
    if (format != 'default') {
      URL <- paste0(URL, "&f=", format)
    }
    session$snapshot <- if (isTRUE(snapshot) 
                           || tolower(substr(snapshot, 1, 1)) == "y") {
      URL <- paste0(URL, "&s=y")
      TRUE
    } else FALSE
    
    #session$call <- match.call() 
    session$URL <- URL
    session$currency.pairs <- currency.pairs
    session$username <- username
    session$password <- password
    session$qualifier <- qualifier
    session$format <- format
    
    session$id <- readLines(URL) #returns the session id
    
    session$last.used <- with(session, 
      if (grepl(paste(username, password, sep=":"), id)) {
        session$active <- TRUE
        Sys.time()
      } else {
        NA
      })
    
    #session$isActive <- function() {
    #  !is.na(id) && !isTRUE(snapshot) && 
    #    difftime(Sys.time(), connected.at, units='secs') <= 60
    #}
    #environment(session$IsOpen) <- as.environment(session)
        
    class(session) <- c("TFXsession", "environment")
    session
  }
}

is.TFXsession <- function(x) {
  inherits(x, 'TFXsession')
}


#' Disconnect a session
#'
#' Disconnect a session (make it inactive).
#'
#' @param x an object to disconnect
#' @param ... other arguments for methods
#' @seealso \code{\link{ConnectTrueFX}}, \code{\link{Reconnect}}
#' @examples
#' \dontrun{
#' sess <- ConnectTrueFX(username='JSTrader', password='Ou812')
#' isActive(sess) #TRUE
#' Disconnect(sess)
#' isActive(sess) #FALSE
#' }
#' @rdname Disconnect
#' @export
Disconnect <- function(x, ...) { UseMethod("Disconnect") }

#' @return a disconnected/inactive \code{TFXsession} object
#'
#' @rdname Disconnect
#' @method Disconnect TFXsession
#' @S3method Disconnect TFXsession
Disconnect.TFXsession <- function(x, ...) {
  stopifnot(inherits(x, 'TFXsession'))
  x$last.used <- NA
  x$active <- FALSE
  readLines(paste0("http://webrates.truefx.com/rates/connect.html?di=", 
                   x$id))
  x
}


#' Is a session active?
#' 
#' Test to see if a session is active
#' 
#' In order for a TFXsession to be active, it must have been authenticated less
#' than 1 minute ago.  If it was created with \code{snapshot=TRUE} it will
#' become inactive after it is used once.
#' @param x an object to test
#' @param ... other arguments for methods
#' @note This function assumes that if the session has not been used in 60 
#' seconds is not active, even though TrueFX(tm) sessions actually stay active 
#' for roughly 70 seconds.
#' @examples
#' \dontrun{
#' sess <- ConnectTrueFX("GBP/JPY", username='JSTrader', password='Ou812')
#' isActive(sess) #TRUE
#' }
#' @rdname isActive
#' @export isActive
isActive <- function(x, ...) { UseMethod("isActive") }

#' @return logical
#' 
#' @rdname isActive
#' @method isActive TFXsession
#' @S3method isActive TFXsession
isActive.TFXsession <- function(x, ...) {
  stopifnot(inherits(x, "TFXsession"))
  # A session will terminate immediately after being used if snapshot == TRUE.
  # A session will also terminate after roughly 1 minute of inactivity. 
  # (actually 70 seconds -- I think)

  isTRUE(x$active) && !is.na(x$last.used) && 
    difftime(Sys.time(), x$last.used, units='secs') <= 60
}


#' Reconnect a session that is no longer active
#' 
#' \code{Reconnect} will create a new session and update the `id` to the new
#' authenticated id returned by the TrueFX(tm) server.
#'
#' After roughly 70 seconds, an authenticated TrueFX(tm) session will time-out.
#' Also, a connection made with \code{snapshot=FALSE} will be disconnected after
#' it is used once.
#'
#' A non-active TrueFX(tm) session id is treated like an unauthenticated 
#' session.
#'
#' @param x an object to be re-connected
#' @param ... other args for methods
#' @examples
#' ## Cannot run because there may not be an internet connection
#' \dontrun{
#' ## You must use your username and password instead of JSTrader and Ou812
#' sess <- ConnectTrueFX("USD/JPY", username='JSTrader', password='Ou812')
#' Disconnect(sess)
#' isActive(sess) #FALSE
#' Reconnect(sess)
#' isActive(sess) #TRUE
#' }
#' @rdname Reconnect
#' @export Reconnect
Reconnect <- function(x, ...) { UseMethod("Reconnect") }

#' @return a \code{TFXsession} object of an active/authenticated session.
#' 
#' @rdname Reconnect
#' @method Reconnect TFXsession
#' @S3method Reconnect TFXsession
Reconnect.TFXsession <- function(x, ...) {
  stopifnot(inherits(x, 'TFXsession'))
  x$id <- readLines(x$URL) #returns the session id
  x$last.used <- if (grepl(paste(x$username, x$password, sep=":"), x$id)) { 
    x$active <- TRUE
    Sys.time() 
  } else { NA }
  x
}


#' Query TrueFX(tm)
#' 
#' Create a session with TrueFX(tm) and request market data.
#' 
#' If no \code{currency.pairs} are provided to \code{ConnectTrueFX}, the 15 
#' pairs for which TrueFX(tm) offers historical data will be used.  Note that 
#' only the first 10 of these are returned in an unauthenticated session.
#' 
#' \code{ConnectTrueFX} will create a \code{TFXsession} classed object that can 
#' be used in calls to \code{QueryTrueFX} to request market data.  
#' 
#' Of the three \code{format}s, \dQuote{default} is the most timely (updates 
#' first)and \dQuote{csv} is the most delayed (updates last)
#' 
#' the \dQuote{csv} and \dQuote{html} formats have the \dQuote{High} and 
#' \dQuote{Low} columns backwards. (\dQuote{default} does not).  This may be 
#' corrected for in a future release if the TrueFX(tm) Web service doesn't 
#' correct it first.
#'
#' @param currency.pairs character vector, or comma delimited string of Symbols
#'   (ISO names) of currency pairs.  (e.g. \code{"EUR/USD,AUD/USD"}, or 
#'   \code{c("EUR/USD", "AUD/USD")}).  If \code{missing} or if 
#'   \code{nchar(currency.pairs) < 1}, the Symbols of all currency pairs for
#'   which TrueFX(tm) provides historical data will be used (see references 
#'   section).
#' @param username character.  A registered TrueFX(tm) user name; required to 
#'   establish an authenticated session.
#' @param password character. A registered TrueFX(tm) password; required to 
#'   establish an authenticated session
#' @param qualifier any string; required to establish an authenticated session.
#'   (\dQuote{default} by default)
#' @param format One of \dQuote{default}, \dQuote{csv}, or \dQuote{html}. 
#'   Indicates the format for the HTTP Response.
#' @param snapshot logical.  No incremental updates if \code{TRUE}
#' 
#' @param session a \code{TFXsession} object created by \code{ConnectTrueFX}.  
#' @param parse.response logical. Should the results be passed through 
#'   \code{ParseTrueFX} before returning?
#' @param pretty logical.  Passed to \code{ParseTrueFX}.  Indicates whether to
#'   format the parsed results and convert to \code{data.frame}. 
#'   Ignored if \code{parse.response} is not \code{TRUE}
#' @param reconnect logical.  If the TFXsession has timed out, should it be
#'   reconnected?
#' @return \code{ConnectTrueFX} returns a \code{TFXsession} object that is a 
#'   TrueFX(tm) server-generated session ID returned with a successful 
#'   authenticated session request.  It is a colon delimited string with 
#'   username, password, qualifier, and the time (in milliseconds) that the 
#'   session was created.  
#'
#'   \code{QueryTrueFX} returns the results of a TrueFX(tm) request using 
#"   \code{TFXsession} object returned by \code{ConnectTrueFX}
#' @author Garrett See
#' @references 
#' \url{http://www.truefx.com/dev/data/TrueFX_MarketDataWebAPI_DeveloperGuide.pdf}
#'
#' \url{http://truefx.com/?page=downloads} to see for which pairs TrueFX(tm) 
#'   provides historical data.
#' @seealso \code{\link{ParseTrueFX}}, \code{\link{Reconnect}}, 
#'   \code{\link{TrueFXRef}}
#' @note the formal arguments start with the same lowercase letter as their 
#'   corresponding TrueFX(tm) Market Data Web Query Parameters
#' @examples
#' ## Cannot run these because there may not be an internet connection
#' \dontrun{
#' QueryTrueFX()  #unauthenticated
#' QueryTrueFX(pretty=FALSE)
#' QueryTrueFX(parse=FALSE)
#' 
#' ## For authenticated session, you must have a username and password (it's free).
#' ## Use your username and passward instead of JSTrader and Ou812
#' id <- ConnectTrueFX('EUR/USD,GBP/USD', username='JSTrader', password='Ou812')
#' QueryTrueFX(id)
#' QueryTrueFX(ConnectTrueFX(username='JSTrader', password='Ou812', 
#'                           format='csv'), parse=FALSE)
#' 
#' QueryTrueFX(ConnectTrueFX(username='JSTrader', password='Ou812', 
#'                           format='html'), parse=FALSE)
#'
#' ## If you have shiny installed 
#' ## install.packages("shiny", repos="http://rstudio.org/_packages")
#' library(shiny)
#' runGist("4122626") 
#' ## view the code for this shiny app at 
#' #browseURL("https://gist.github.com/4122626")
#' }
#' @export
#' @rdname QueryTrueFX
QueryTrueFX <- function(session, parse.response=TRUE, pretty=TRUE, 
    reconnect=TRUE) {
  if (missing(session)) {
    if (isTRUE(parse.response)) {
      return(ParseTrueFX(readLines(
        "http://webrates.truefx.com/rates/connect.html"), pretty=pretty))
    } else return(readLines("http://webrates.truefx.com/rates/connect.html")) 
  }
  if (!inherits(session, "TFXsession")) {
    stop("session is not a TFXsession object created by ConnectTrueFX")
    # or should it warn and
    # return(readLines("http://webrates.truefx.com/rates/connect.html"))
  }
  if (session$id == "not authorized") stop("not authorized")
  if (!isActive(session)) {
    if (isTRUE(reconnect)) {
      #warning("session is no longer active. Reconnecting ...")
      session <- Reconnect(session)
    } else stop("'session' is not connected and 'reconnect' is not TRUE")
  } 
  if (isTRUE(session$snapshot)) {
    session$active <- FALSE
  }
  session$last.used <- Sys.time() 
  ## request
  if (isTRUE(parse.response)) {
    return(ParseTrueFX(readLines(paste0(
      "http://webrates.truefx.com/rates/connect.html?id=", session$id)), 
                       pretty=pretty))
  }
  readLines(paste0("http://webrates.truefx.com/rates/connect.html?id=", 
                   session$id))
  # The next line would disconnect
  #readLines(paste0("http://webrates.truefx.com/rates/connect.html?di=", 
  #                 session$id))
}

#' Parse TrueFX(tm) response
#' 
#' Parse the results of a TrueFX(tm) query.  
#' 
#' This function will parse the results of a call to \code{\link{QueryTrueFX}}. 
#' It can handle any of the three TrueFX(tm) response formats: \dQuote{default}, 
#' \dQuote{csv}, or \dQuote{html}.  By default, it will convert the results 
#' into a nicely formatted \code{data.frame}.  If, called with 
#' \code{pretty=FALSE}, a list of strings will be returned.
#' 
#' All times are in GMT
#' 
#' @param x The response from a TrueFX(tm) request.  Can be any of the three
#'   formats: \code{default}, \code{csv} or \code{html}
#' @param pretty logical. If \code{TRUE} (Default), output will be converted to 
#'   a data.frame and columns will be converted from character to the 
#'   appropriate classes and combined.
#' @return By default, a \code{data.frame} is returned that has columns 
#'   \dQuote{Bid.Price}, \dQuote{Ask.Price}, \dQuote{High}, \dQuote{Low},
#'   and \dQuote{TimeStamp}.  If called with \code{pretty=FALSE}, a list of 
#'   character vectors -- named \dQuote{Symbol}, \dQuote{BidBigNumber}, 
#'   \dQuote{BidPip}, \dQuote{OfferBigNumber}, \dQuote{OfferPip}, 
#'   \dQuote{High}, \dQuote{Low}, \dQuote{TimeStamp} -- will be returned.
#'   
#'   If the format is \dQuote{html}, there will also be an \dQuote{Open} column
#' @author Garrett See
#' @references 
#' \url{http://www.truefx.com/dev/data/TrueFX_MarketDataWebAPI_DeveloperGuide.pdf}
#' @seealso \code{\link{QueryTrueFX}}, \code{\link{TrueFXRef}}
#' @note Although the TrueFX(tm) Market Data Web API Developer Guide indicates 
#'   that both the \dQuote{csv} and \dQuote{html} formats include values for 
#'   \dQuote{Open}, only the \dQuote{html} format actually does.
#' @examples 
#' # x <- QueryTrueFX()  #Cannot run this if no internet connection
#' x <- paste0("EUR/USDUSD/JPY1.31#81.9085661.31#81.9435941.31990#81.6421.3182",
#'             "1#81.50413351311514701335131150004")
#' ParseTrueFX(x)
#' ParseTrueFX(x, pretty=FALSE)
#' @importFrom XML readHTMLTable
#' @export
ParseTrueFX <- function(x, pretty=TRUE) {
  PasteFigurePip <- function(figure, pip) {
    out <- gsub("\\.", "", paste0(sprintf("%04s", as.numeric(figure)), 
                                  sprintf("%03s", as.numeric(pip))))
      tmp <- sprintf("%04s", as.numeric(figure))
      loc <- -grep("\\.", tmp)
      # if it doesn't have a dot, add one at the end
      tmp[loc] <- paste0(tmp[loc], ".")
    as.numeric(paste0(tmp, gsub(" ", 0, sprintf("%03s", as.numeric(pip)))))
  }
  
  if (any(grepl(",", x))) {  # It's in csv format
    x <- x[x != ""]
    if (!isTRUE(pretty)) {
      return(as.list(read.csv(text=x, header=FALSE, stringsAsFactors=FALSE, 
                              col.names = c("Symbol", "TimeStamp", 
                                            "BidBigNumber", 
                                            "BidPip", "OfferBigNumber", 
                                            "OfferPip", "High", "Low", "Open"),
                              colClasses = 'character')))
    } else {
      tmp <- read.csv(text=x, header=FALSE, stringsAsFactors=FALSE, 
                      col.names = c("Symbol", "TimeStamp", "BidBigNumber", 
                                    "BidPip", "OfferBigNumber", "OfferPip", 
                                    "High", "Low", "Open"))
      return(data.frame(Symbol = tmp[["Symbol"]],
        Bid.Price = PasteFigurePip(tmp[["BidBigNumber"]], tmp[["BidPip"]]),
        Ask.Price = PasteFigurePip(tmp[["OfferBigNumber"]], tmp[["OfferPip"]]),
        High = tmp[["High"]],
        Low = tmp[["Low"]],
        Open = tmp[["Open"]],
        TimeStamp = as.POSIXct(as.numeric(tmp[["TimeStamp"]]) / 1000, 
                               origin='1970-01-01', tz='GMT'), 
        stringsAsFactors=FALSE))
    }
  } else if (substr(x, 1, 7) == "<table>") { # It's an HTML table
    out <- readHTMLTable(x, as.data.frame=FALSE)[[1]]
    names(out) <- c("Symbol", "TimeStamp", "BidBigNumber", "BidPip", 
                    "OfferBigNumber", "OfferPip", "High", "Low", "Open")
    if (!isTRUE(pretty)) {
      return(out)
    } else {
      return(data.frame(Symbol=out[["Symbol"]],
                        Bid.Price=as.numeric(paste0(out[["BidBigNumber"]],
                                                    out[["BidPip"]])),
                        Ask.Price=as.numeric(paste0(out[["OfferBigNumber"]],
                                                    out[["OfferPip"]])),
                        High=as.numeric(out[["High"]]),
                        Low=as.numeric(out[["Low"]]),
                        Open=as.numeric(out[["Open"]]),
                        TimeStamp=as.POSIXct(as.numeric(out$TimeStamp) / 1000, 
                                             origin='1970-01-01', tz='GMT'),
                        stringsAsFactors=FALSE))
    }
  }
  # Otherwise, it's a concatenated string
  npairs <- nchar(gsub("[0-9.#]", "", x)) / 7
  .ReadSection <- function(string, by) {
    end <- by * npairs
    if (end > 0) substring(string, seq(1, end, by), seq(by, end, by))
  }
  # See page 3 of
  #http://www.truefx.com/dev/data/TrueFX_MarketDataWebAPI_DeveloperGuide.pdf
  Ns <- c(7, 4, 3, 4, 3, 7, 7, 13)
  ep <- c(0, cumsum(npairs * Ns))
  out <- lapply(1:(length(ep) - 1), function(i) {
    beg <- (ep[i] + 1)
    end <- ep[i + 1]
    .string <- substr(x, beg, end)
    .ReadSection(.string, Ns[i])
  })
  names(out) <- c("Symbol", "BidBigNumber", "BidPip", "OfferBigNumber", 
                  "OfferPip", "High", "Low", "TimeStamp")
  if (!isTRUE(pretty)) {
    return(out)
  }
  out <- lapply(out, gsub, pattern="#", replacement=0)  
  data.frame(Symbol=out[["Symbol"]],
             Bid.Price=as.numeric(paste0(out[["BidBigNumber"]],
                                         out[["BidPip"]])),
             Ask.Price=as.numeric(paste0(out[["OfferBigNumber"]],
                                         out[["OfferPip"]])),
             High=as.numeric(out[["High"]]),
             Low=as.numeric(out[["Low"]]),
             TimeStamp=as.POSIXct(as.numeric(out$TimeStamp) / 1000, 
                                  origin='1970-01-01', tz='GMT'),
             stringsAsFactors=FALSE)
}


