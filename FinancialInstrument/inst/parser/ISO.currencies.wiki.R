#' Define currencies using the tables found on oanda's website.
#'
#' If you do not provide \code{Symbols} all oanda curriencies will be defined.
#' If you do provide \code{Symbols} only the Symbols you provided will be defined. 
#' @param Symbols
#' @param silent
#' @return the names of the currecies that were defined.  Called for side-effect
#' @references \url{http://www.oanda.com/help/currency-iso-code}
#' @author Garrett See
#' @examples
#' \dontrun{
#' define_currencies.oanda(c("EUR","GBP","JPY"))
#' define_currencies.oanda()
#' }
define_currencies.oanda <- function(Symbols, silent=FALSE) {
    if (!("package:XML" %in% search() || require("XML",quietly=TRUE))) 
        stop("Please install the XML package before using this function.")
    x <- readHTMLTable("http://www.oanda.com/help/currency-iso-code")
    x <- lapply(x, function(xx) xx[-1,])
    #all.syms <- unname(do.call(c,lapply(x, function(xx) as.character(xx[,1]))))
    df <- do.call(rbind, lapply(x, function(xx) cbind(as.character(xx[,1]), as.character(xx[,2]))))
    if(missing(Symbols)) Symbols <- df[,1]
    df <- df[df[,1] %in% Symbols,]
    apply(df,1,function(X) currency(X[1], description=X[2],defined.by='oanda'))
}


#http://en.wikipedia.org/wiki/ISO_4217#Active_codes

#' Define currency instruments using the tables found on the ISO_4217 wikipedia page
#'
#' If you do not provide \code{Symbols} all active ISO 4127 curriencies will be defined.
#' If you do provide \code{Symbols} only the Symbols you provided will be defined. Also,
#' if you provide some \code{Symbols} that are not active ISO currencies, it will try to 
#' find them in the "Without currency code" and "Historic currency codes" tables.
#' @param Symbols
#' @param silent
#' @return the names of the currecies that were defined.  
#' If \code{Symbols} was provided, they will be in the order they were found (ISO, non-ISO, historic). 
#' Called for side-effect
#' @references \url{http://en.wikipedia.org/wiki/ISO_4217}
#' @examples
#' \dontrun{
#' define_currencies.wiki(c("USD","EUR","ADP","ETB","GBP","BTC"))
#' define_currencies.wiki()
#' }
define_currencies.wiki <- function(Symbols, silent=FALSE) {
    if (!("package:XML" %in% search() || require("XML",quietly=TRUE))) 
        stop("Please install the XML package before using this function.")
    x <- readHTMLTable("http://en.wikipedia.org/wiki/ISO_4217")
    #"The following is a list of active codes of official ISO 4217 currency names."        
    ccy <- x[[2]] #active ISO 4217 currencies
    if (!missing(Symbols)) {
        ccy <- rbind(ccy, x[[3]]) #add non-ISO... things like Bitcoin
        ccy <- ccy[ccy$Code %in% Symbols,]
    } else Symbols <- NULL
    out <- unname(apply(ccy, 1, function(xx) currency(xx[1], identifiers=list(Num=xx[2]), digits.after.dec=xx[3], 
                                                        description=xx[4], country=xx[5], defined.by='wiki')))
    if (!is.null(Symbols) && !identical(character(0), Symbols[!Symbols %in% ccy$Code])) {
        if (!silent) warning(paste("The following are historical,",
                "and are no longer active:", Symbols[!Symbols %in% ccy$Code]))
        hccy <- x[[4]]
        hccy <- hccy[hccy$Code %in% Symbols,]
        out <- c(out, unname(apply(hccy, 1, function(xx) {
                    currency(xx[1], description=xx[4], used.from=xx[5], used.until=xx[6], replaced.by=xx[7], defined.by='wiki')
                 })))
    }
    out
}

#rm_currencies()
#define_currencies.wiki(c("USD","XAU"))
#exchange_rate("XAUUSD",src=list(src='oanda',name='XAU/USD'))

#define_currencies.wiki(c("JPY","EUR","ADP","ETB","GBP","BTC"))
#define_currencies.wiki()

#_______________________________________________________________________________________

#x <- readHTMLTable("http://en.wikipedia.org/wiki/ISO_4217")

#"A number of territories are not included in ISO 4217, because their currencies are: 
#(a) not per se an independent currency but a variant of another currency, 
#(b) a legal tender only issued as commemorative banknotes or coinage, or 
#(c) a currency of an unrecognized or partially recognized state. These currencies are:"
#nonISO.ccy <- x[[3]]

#"A number of currencies were official ISO 4217 currency codes and currency names 
#until their replacement by the euro or other currencies. 
#The table below shows the ISO currency codes of former currencies and their common names 
#(which do not always match the ISO 4217 name)."
#historic.ccy <- x[[4]]

