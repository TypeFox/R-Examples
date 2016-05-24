###############################################################################
# R (http://r-project.org/) Instrument Class Model
#
# Copyright (c) 2009-2012
# Peter Carl, Dirk Eddelbuettel, Jeffrey Ryan, 
# Joshua Ulrich, Brian G. Peterson, and Garrett See
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: find.instrument.R 1092 2012-06-30 00:06:09Z gsee $
#
###############################################################################


#' Find the primary_ids of instruments that contain certain strings
#' 
#' Uses regular expression matching to find \code{\link{instrument}}s
#' 
#' @param text character string containing a regular expression.  This is used 
#'   by \code{\link{grep}} (see also) as the \code{pattern} argument.
#' @param where if \dQuote{anywhere} all levels/attributes of the instruments
#' will be searched.  Otherwise, \code{where} can be used to specify in which
#' levels/attributes to look. (e.g. \code{c("name", "description")} would only
#' look for \code{text} in those 2 places.
#' @param Symbols the character ids of instruments to be searched. All are 
#' are searched by default.
#' @param ignore.case passed to \code{\link{grep}}; if \code{FALSE}, the pattern
#' matching is case sensitive and if \code{TRUE}, case is ignored during 
#' matching.
#' @param exclude character vector of names of levels/attributes that should not
#' be searched.
#' @param ... other arguments to pass through to \code{\link{grep}}
#' @return character vector of primary_ids of instruments that contain the 
#' sought after \code{text}.
#' @author Garrett See
#' @seealso \code{\link{buildHierarchy}}, \code{\link{instrument.table}}, 
#' \code{\link{regex}}
#' @examples
#' \dontrun{
#' instruments.bak <- as.list(FinancialInstrument:::.instrument, all.names=TRUE)
#' rm_instruments(keep.currencies=FALSE)
#' currency("USD")
#' stock("SPY", "USD", description="S&P 500 ETF")
#' stock("DIA", "USD", description="DJIA ETF")
#' stock(c("AA", "AXP", "BA", "BAC", "CAT"), "USD", members.of='DJIA')
#' stock("BMW", currency("EUR"))
#' find.instrument("ETF")
#' find.instrument("DJIA") 
#' find.instrument("DJIA", "members.of")
#' find.instrument("USD")
#' find.instrument("EUR")
#' find.instrument("EUR", Symbols=ls_stocks())
#' find.instrument("USD", "type")
#'
#' ## Can be combined with buildHierachy
#' buildHierarchy(find.instrument("ETF"), "type", "description")
#'
#' ## Cleanup. restore previous instrument environment
#' rm_instruments(); rm_currencies()
#' loadInstruments(instruments.bak)
#' }
#' @export
find.instrument <- function(text, where='anywhere', Symbols = ls_instruments(),
                            ignore.case=TRUE, exclude=NULL, ...) {
    tbl <- if (length(where) == 1 && where == "anywhere") {
        instrument.table(Symbols, exclude=exclude)
    } else buildHierarchy(Symbols, where[!where %in% exclude])
    unique(tbl[unique(unname(unlist(apply(tbl, 2, function(x) 
        grep(pattern=text, x=x, ignore.case=ignore.case, 
             useBytes=TRUE, ...))))), 1])
}
