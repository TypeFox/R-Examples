#' update iShares and SPDR ETF metadata
#' 
#' This will update previously defined iShares or SPDR ETF \code{instrument}s. 
#' Both functions will add attributes for \dQuote{Name}, and \dQuote{FundFamily}
#' (\dQuote{iShares} or \dQuote{SPDR}). \code{update_instruments.iShares} will 
#' also add an attribute for \dQuote{MgmtFees}
#' 
#' @param Symbols character vector of iShares ETF ticker symbols.  If not 
#'   specified, \code{unique(c(ls_funds(), ls_stocks()))} will be used.
#' @param silent silence the warning that no iShares are defined? 
#' @return called for side-effect
#' @author Garrett See
#' @seealso \code{update_instruments.yahoo}, \code{update_instruments.TTR},
#'   \code{twsInstrument:::update_instruments.IB}, 
#'   \code{update_instruments.instrument}, 
#'   \code{\link{update_instruments.morningstar}},
#'   \code{\link{update_instruments.masterDATA}}
#' @note \code{update_instruments.SPDR} will probably NOT work on Windows 
#'   because in the call to \code{download.file} it uses \code{method=curl} 
#'   since it has to download from an https URL scheme.
#' @references \url{http://us.ishares.com/home.htm}, 
#'   \url{https://www.spdrs.com/}
#' @examples
#' \dontrun{
#' stock("IWC", currency("USD"))
#' update_instruments.iShares("IWC")
#' getInstrument("IWC")
#' 
#' Symbols <- stock(c("SPY", "JNK"), currency("USD"))
#' update_instruments.SPDR(Symbols)
#' buildHierarchy(c("SPY", "JNK"), "Name")
#' }
#' @export
#' @rdname update_instruments.iShares
update_instruments.iShares <- function(Symbols, silent=FALSE) {
  tmp <- tempfile()
  lnk <- paste("http://us.ishares.com/product_info/fund/excel_returns.htm",
                "?assetClassCd=EQ&ticker=&asofDt=", sep="")
  download.file(lnk, destfile=tmp, quiet=TRUE)
  fr <- read.csv(tmp, skip=3, stringsAsFactors=FALSE, header=FALSE)
  colnames(fr) <- read.delim(text=readLines(tmp, 1), sep=",", header=FALSE, 
                             stringsAsFactors=FALSE)
  unlink(tmp)  
  tickers <- fr[["Ticker"]]
  if (missing(Symbols)) {
    Symbols <- unique(c(ls_funds(), ls_stocks()))
  }
  s <- Symbols[Symbols %in% tickers]
  if (length(s) == 0) {
    if (!isTRUE(silent)) {
      warning('iShares must be defined before this can update them.')
    }
    return(invisible())
  }
  dat <- fr[fr$Ticker %in% s, ]
  for(i in seq_len(NROW(dat))) {
    instrument_attr(dat[i, "Ticker"], "MgmtFees", dat[i, "Mgmt Fees"])
    instrument_attr(dat[i, "Ticker"], "Name", dat[i, "iShares Fund Name"])
    instrument_attr(dat[i, "Ticker"], "FundFamily", "iShares")
  } 
  invisible(lapply(s, function(x) {
    db <- getInstrument(x)[["defined.by"]]
    db <- if (is.null(db)) {
      "iShares"
    } else paste(unique(c(strsplit(db, ";")[[1]], 
              "iShares")), collapse=";")
    instrument_attr(x, "defined.by", db)
    instrument_attr(x, "updated", Sys.time())
  }))
  return(s)
}


#' @export 
#' @rdname update_instruments.iShares
update_instruments.SPDR <- function(Symbols, silent=FALSE) {
  tmp <- tempfile()
  lnk <- paste("https://www.spdrs.com/library-content/public/public-files/etf",
                "nav.csv?docname=Most+Recent+Net+Asset+Values&onyx_code1=1299",
               sep="")
  download.file(lnk, destfile=tmp, method="curl")
  fr <- read.csv(tmp, skip=1, stringsAsFactors=FALSE)
  DATE <- gsub(" ", "", sub("DATE,", "", readLines(tmp, 1)))
  unlink(tmp)
  tickers <- fr[["TICKER"]] <- gsub("^\\s+|\\s+$", "", fr[["TICKER"]])
  if (missing(Symbols)) {
    Symbols <- unique(c(ls_funds(), ls_stocks()))
  }
  s <- Symbols[Symbols %in% tickers]
  if (length(s) == 0) {
    if (!isTRUE(silent)) {
      warning('instruments must be defined before this can update them.')
    }
    return(invisible())
  }
  dat <- fr[fr[["TICKER"]] %in% s, ]
  for (i in seq_len(NROW(dat))) {
    instrument_attr(dat[i, "TICKER"], "Name", 
                    gsub("^\\s+|\\s+$", "", dat[i, "NAME"]))
    instrument_attr(dat[i, "TICKER"], "FundFamily", "SPDR")
    add.identifier(dat[i, "TICKER"], 
                   SPDR=gsub("^\\s+|\\s+$", "", dat[i, "FUND"]))
    add.identifier(dat[i, "TICKER"], CUSIP=gsub("\'", "", fr[i, "CUSIP"]))
  }
  invisible(lapply(s, function(x) {
    db <- getInstrument(x)[["defined.by"]]
    db <- if (is.null(db)) {
      "SPDR"
    } else paste(unique(c(strsplit(db, ";")[[1]], "SPDR")), collapse=";")
    instrument_attr(x, "defined.by", db)
    instrument_attr(x, "updated", Sys.time())
  }))
  return(s)
}

