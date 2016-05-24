#' Update instrument metadata for ETFs
#' 
#' Currently, this only updates ETFs.  It will add \dQuote{msName} and
#' \dQuote{msCategory} attributes to the instruments. (ms for morningstar)
#' @param Symbols character vector of Symbols of ETFs
#' @param silent silence warnings?
#' @return called for side-effect. 
#' @author Garrett See
#' @references \url{http://www.morningstar.com}
#' @seealso \code{\link{update_instruments.yahoo}}, 
#' \code{\link{update_instruments.TTR}}
#' \code{\link{update_instruments.iShares}}
#' @examples
#' \dontrun{
#' ## backup .instrument environment
#' ibak <- as.list(FinancialInstrument:::.instrument) 
#' rm_instruments()
#' stock(s <- c("SPY", "USO", "LQD"), currency("USD"))
#' update_instruments.morningstar(s)
#' instrument.table(s)
#' ## cleanup and restore instrument environment
#' rm_instruments(keep.currencies=FALSE)
#' loadInstruments(ibak)
#' }
#' @export
update_instruments.morningstar <- function(Symbols, silent=FALSE) {
    require(XML)
    x <- readHTMLTable(paste("http://news.morningstar.com/etf/Lists/ETFReturn",
                             "s.html?topNum=All&lastRecNum=1000&curField=8&ca",
                             "tegory=0", sep=""), stringsAsFactors=FALSE)
    x <- x[[which.max(sapply(x, NROW))]]
    colnames(x) <- x[1, ]
    x <- x[-c(1:2), -1]
    x <- x[!is.na(x[, 1]), ]
    x <- x[!duplicated(x[, 1]), ]
    tickers <- gsub(".*\\(|*\\)", "", x[,1])
    x <- x[tickers != "", ]
    tickers <- tickers[tickers != ""]
    rownames(x) <- tickers
    if (missing(Symbols)) {
        Symbols <- unique(c(ls_funds(), ls_stocks()))
    }
    s <- Symbols[Symbols %in% tickers]
    if (length(s) > 0) {
        # only those that inherit stock or fund
        s <- s[sapply(lapply(s, getInstrument, type=c("stock", "fund"), 
                             silent = TRUE), is.instrument)]
    }
    if (length(s) == 0) {
        if (!isTRUE(silent)) {
            warning("instruments must be defined before this can update them.")
        }
        return(invisible())
    }
    x <- x[rownames(x) %in% s, ]
    rn <- rownames(x)
    for (i in 1:NROW(x)) {
        instrument_attr(rn[i], "msName", x$Name[i])
        instrument_attr(rn[i], "msCategory", x$Category[i])
        #instrument_attr(x$Symbol[i], "msTradingVolume", 
        #                as.numeric(gsub(",", "", x$TradingVolume[i])))
        db <- getInstrument(rn[i])$defined.by
        instrument_attr(rn[i], "defined.by", paste(c(db, "morningstar"), 
                                                   collapse=";"))
        instrument_attr(rn[i], "updated", Sys.time())
    }
    return(s)
}


#' @export
#' @rdname update_instruments.morningstar
update_instruments.ms <- update_instruments.morningstar

