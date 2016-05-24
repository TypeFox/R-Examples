#' get the components of the Dow Jones Industrial Average
#'
#' download a data.frame of the 30 Dow Jones components.
#' @return 30 by 5 data.frame with columns \sQuote{Symbol}, \sQuote{Name}, \sQuote{Last.Trade}, \sQuote{Change}, \sQuote{Volume}
#' @references \url{'http://finance.yahoo.com/q/cp?s=^DJI+Components'}
DJIcomponents <- function() {
    if (!("package:XML" %in% search() || require("XML",quietly=TRUE))) {
        stop("Please install the XML package before using this function.")
    }
    djicomp <- readHTMLTable('http://finance.yahoo.com/q/cp?s=^DJI+Components')
    data.frame(djicomp[[grep("Symbol", djicomp)]])
}

#' fetch the current divisor for the Dow Jones Industrial Average from Barrons
#' @return numeric
#' @references \url{'http://online.barrons.com/mdc/public/page/9_3022-djiahourly.html?mod=mdc_h_usshl'}
dow.divisor <- function() {
    wp <- readLines('http://online.barrons.com/mdc/public/page/9_3022-djiahourly.html?mod=mdc_h_usshl')
    #wp2 <- wp[grep("30 INDUSTRIALS:",wp)]
    #as.numeric(gsub(')</td>','',strsplit(wp2, 'divisor: ')[[1]][2]))
    as.numeric(gsub(')</td>','',strsplit(wp[grep("30 INDUSTRIALS:",wp)], 'divisor: ')[[1]][2]))
}

currency('USD')
DJIA.members <- stock(DJIcomponents()$Symbol, currency="USD", member.of='DJIA')
getSymbols(DJIA.members) # <-- getting data first is not required, but may be preferable
synthetic.instrument("DJIA","USD", members=DJIA.members, memberratio=rep(1,length(DJIA.members)), 
                     multiplier=1/dow.divisor(), tick_size=0.01, description='Dow Jones Industrial Average')
buildBasket('DJIA') #theoretical index
tail(DJIA)

getSymbols("^DJI") #actual index
tail(DJI)

