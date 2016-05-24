
#' Define an Index and it's components using yahoo
#'
#' Get the components of an index and define instruments
#' @param Symbol character yahoo ticker symbol for a stock index (e.g. "^DJI" or "^GDAXI")
#' @param currency
#' @return called for side-effect, but it will return a list with 2 components:
#' \item{synthetic}{name of the \code{\link{synthetic}} instrument that was defined to hold the metadata of the index.} 
#' \item{stock}{name of the component \code{\link{stock}}s that were defined}
#' @note Depends on XML package
#' @author Garrett See
#' @examples
#' \dontrun{
#' define_components.yahoo('^STOXX50E', 'EUR')
#' define_components.yahoo('^DJI', 'USD')
#' }
define_components.yahoo <- function(Symbol, currency) {
    require(FinancialInstrument)
    require(XML)
    ccy <- currency(currency) #make sure it's defined
    x <- readHTMLTable(paste("http://finance.yahoo.com/q/cp?s=",Symbol,"+Components", sep=""))
    mdata <- x[which.max(sapply(x, NROW))][[1]][,1:2]
    new.Symbol <- synthetic(Symbol, ccy, members=paste(mdata[,1]), src=list(src='yahoo',name=Symbol), identifiers=list(yahoo=Symbol))
    if (!identical(integer(0), grep("There is no Components data", mdata[,1]))) stop("No Components Data Available for ", Symbol)
    stks <- rep(NA_character_, NROW(mdata))    
    for (i in 1:NROW(mdata)) {
        tmpsym <- paste(mdata[i,1])        
        stks[i] <- stock(tmpsym, currency=ccy, Name=paste(mdata[i,2]), member.of=new.Symbol)
        if (!identical(tmpsym,stks[i])) { 
            # add info about how to find data and metadata
            instrument_attr(stks[i], 'src', list(src='yahoo', name=tmpsym))           
            instrument_attr(stks[i], 'identifiers', list(yahoo=tmpsym)) 
        } 
    }
    list(synthetic=new.Symbol, stock=stks)
}

#define_components.yahoo('^STOXX50E','EUR')
#Symbol <- '^STOXX50E'
#instr <- getInstrument("^STOXX50E")
#instr
#memb5 <- getInstrument(instr$members[5])
#memb5
#getInstrument(memb5$member.of)


#define_components.yahoo("^GDAXI","EUR")
#getSymbols("X63DU.DE")
#getSymbols("63DU.DE")
#getInstrument("X63DU.DE")
#getInstrument("63DU.DE")




