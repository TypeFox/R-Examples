# Daily historical futures data for 40 products since 1995.

library(FinancialInstrument)

#' get historical futures data from tradingblox.com
#'
#' \code{get_tblox} will download all data available from tradingblox.
#' \code{getSymbols.tblox} will only return the Symbols you ask for.
get_tblox <- function(env='.GlobalEnv') {
# a function to download all data for all 40 instruments
    tmpdir <- tempdir()
    tblox.tmp <- paste(tmpdir, "tblox", sep="/")
    if(!file.exists(tblox.tmp)) dir.create(tblox.tmp)
    tmp <- tempfile()
    download.file("http://www.tradingblox.com/Data/DataOnly.zip", tmp)
    unzip(tmp,exdir=tblox.tmp)
    def <- read.csv('http://www.tradingblox.com/tradingblox/CSIUA/FuturesInfo.txt',skip=1,header=FALSE)
    badsym <- NULL    
    for (i in 1:length(def[,1])){
        if (file.exists(paste(tblox.tmp,'Futures',def[i,4],sep='/'))) {
            dat <-read.csv(paste(tblox.tmp,"/Futures/",def[i,4],sep=""),header=FALSE)
            idx <- as.Date(dat[,1],format='%Y%m%d', origin='1970-01-01')
            x <- xts(dat[,2:9],order.by=as.Date(paste(dat[,1]),format="%Y%m%d"))
            cn <- c(paste('Adj',c("Open","High","Low","Close"),sep="."),'Volume','OpInt','ExpMth','Unadj.Close')
            colnames(x) <- paste(def[i,1],cn,sep=".")        
            assign(paste(def[i,1]), x, pos=env)
        } else badsym <- c(badsym, paste(def[i,1]))
    }
    unlink(tmp)
    unlink(tblox.tmp, recursive=TRUE)
    out <- paste(def[,1])
    out[!out %in% badsym]
}

# a getSymbols method to get only the symbols you specify 
# (it still has to download all the data, but will only read/save data for the symbols you specify)
getSymbols.tblox <-
function (Symbols, env, return.class = "xts", ...) 
{
    #importDefaults("getSymbols.tblox")
    this.env <- environment()
    for (var in names(list(...))) {
        assign(var, list(...)[[var]], this.env)
    }
    default.return.class <- return.class
    if (missing(verbose)) 
        verbose <- FALSE
    if (missing(auto.assign)) 
        auto.assign <- TRUE
    if (!auto.assign) stop("must use auto.assign=TRUE for src='tblox'") 

    tmpdir <- tempdir()
    tblox.tmp <- paste(tmpdir, "tblox", sep="/")
    if(!file.exists(tblox.tmp)) dir.create(tblox.tmp)
    tmp <- tempfile()
    download.file("http://www.tradingblox.com/Data/DataOnly.zip",tmp)
    unzip(tmp,exdir=tblox.tmp)
    def <- read.csv('http://www.tradingblox.com/tradingblox/CSIUA/FuturesInfo.txt',
                    skip=1, header=FALSE, stringsAsFactors=FALSE)
    if (is.null(Symbols) || is.na(Symbols) || Symbols == "all" || Symbols == "")
        Symbols <- def
    sym.out <- NULL
    for (i in match(Symbols, paste(def[,1])) ){
        if (file.exists(paste(tblox.tmp,'Futures',def[i,4],sep='/'))) {
            if (verbose) 
                cat("loading ", def[i, 1], ".....")
            return.class <- getSymbolLookup()[[paste(def[i,1])]]$return.class
            return.class <- ifelse(is.null(return.class), default.return.class, 
                return.class)
            dat <-read.csv(paste(tblox.tmp,"/Futures/",def[i,4],sep=""),header=FALSE)
            fr <- if (verbose) 
                cat("done.\n")
            idx <- as.Date(dat[,1],format='%Y%m%d')
            x <- xts(dat[,2:9],order.by=as.Date(paste(dat[,1]),format="%Y%m%d"))
            cn <- c(paste('Adj',c("Open","High","Low","Close"),sep="."),'Volume','OpInt','ExpMth','Unadj.Close')
            colnames(x) <- paste(def[i,1],cn,sep=".")
            x <- quantmod:::convert.time.series(fr = x, return.class = return.class)
            assign(paste(def[i,1]), x, pos=env)
            sym.out <- c(sym.out, paste(def[i,1]))
        }    
    } 
    unlink(tmp)
    unlink(tblox.tmp, recursive=TRUE)
    return(sym.out)
}

#http://www.tradingblox.com/tradingblox/documentation.htm

#' Define futures with tradingblox data
#'
#' Define futures using the tradingblox.com futures dictionary.
#' @param verbose be verbose?
#' @return called for side-effect
#' @author Garrett See
#' @examples \dontrun{define_futures.tblox()}
define_futures.tblox <- function(verbose=TRUE){
# a function to define metadata for all the futures that are available from tblox
    def <- read.csv('http://www.tradingblox.com/tradingblox/CSIUA/FuturesInfo.txt',skip=1,header=FALSE,stringsAsFactors=FALSE)
    for (i in 1:length(def[,1])) {    
        ccy <- try(getInstrument(paste(def[i,8]),silent=TRUE),silent=TRUE)
        if (inherits(ccy, 'try-error') || !inherits(ccy,'currency')) {
            currency(paste(def[i,8]))
            if (verbose) warning(paste("Created currency", def[i,8], "because it did not exist."))
        }
        tick <- def[i,14]
        if (length(strsplit(tick,"/")[[1]]) == 2) {
            numer <- as.numeric(strsplit(tick,"/")[[1]][1])
            denom <- as.numeric(gsub("h","",strsplit(tick,"/")[[1]][2]))
            tick <- numer / denom
        }
        tmonths <- paste(strsplit(def[i,6],"")[[1]], collapse=",")
        primary_id <- paste(def[i,1])
        instr <- try(getInstrument(primary_id,silent=TRUE),silent=TRUE)
        if (inherits(instr,'try-error') || !is.instrument(instr)) {
            future(primary_id = primary_id,
                    currency = paste(def[i,8]),
                    multiplier = as.numeric(gsub(",","",def[i,9])),
                    tick_size = as.numeric(tick),
                    src = 'tblox',
                    defined.by = 'tblox' )
        }
        instr <- getInstrument(primary_id)
        instr$currency = paste(def[i,8])
        instr$multiplier = as.numeric(gsub(",","",def[i,9]))
        instr$tick_size = as.numeric(tick)
        instr$identifiers = NULL
        instr$description = paste(def[i,2])
        instr$exchange = paste(def[i,5])
       #instr$TradingMonths = tmonths
        instr$month_cycle = tmonths
        instr$ContractSize = paste(def[i,7])
        instr$Margin = as.numeric(gsub(",","",def[i,10]))
        instr$ProductGroup = paste(def[i,17])
        instr$DisplayDigits = paste(def[i,21])
        if (is.null(instr$defined.by)) instr$defined.by <- 'tblox'        
        if(instr$defined.by != 'tblox') 
            instr$defined.by <- paste(c(instr$defined.by, "tblox"), collapse = ";")
        instr$updated <- Sys.time()
        assign(primary_id, instr, pos=FinancialInstrument:::.instrument)
    }
    paste(def[,1])
}

#define_futures.tblox()

#1: symbol
#2: description <- c(description, description)
#5: exchange
#6: TradingMonths
#7: ContractSize
#8: currency
#9: BigPoint Value of a 1.0 movement
#10: Margin
#11: CloseCorrel 
#12: LooseCorrel

#14: tick_unit
#15: minimum_tick
#17: Group
#21: DisplayDigits


#' remove all but adjusted OHLC columns from tblox data
#' @param env environment that holds tblox data
#' @return called for side-effect. 
#' @note must have called define_futures.tblox first
#' @examples
#' \dontrun{
#' get_tblox()
#' define_futures.tblox()
#' reformat_tblox()
#' }
reformat_tblox <- function(symbols, env=.GlobalEnv) {
    f <- function(x){
        nx <- try(get(x,pos=env),silent=TRUE)
        if (!inherits(nx,'try-error')) {
            nx <- nx[,1:5]
            colnames(nx) <- paste(x,c("Open","High","Low","Adjusted.Close","Volume"),sep=".")
            assign(x,nx,pos=env)
        }
    }
    syms <- if(missing(symbols) || is.null(symbols)) {
                if (!require('twsInstrument')) stop('symbols must be supplied if twsInstrument package is not installed.')
                ls_instruments_by("defined.by","tblox")
            } else symbols
    lapply(syms, FUN=f)
    syms
}



