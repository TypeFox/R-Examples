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
# $Id: instrument.table.R 949 2012-02-29 14:31:10Z gsee $
#
###############################################################################

#' Create data.frame with attributes of all instruments
#' 
#' A wrapper for \code{\link{buildHierarchy}}, that defaults to returning all attributes.
#' By default it looks for the instrument with the most attribute levels, and uses those attributes
#' for columns.  If you would prefer to use the attribute levels of a given instrument to build the columns, 
#' use \code{attrs.of}.
#'
#' if there are some attributes that you do not want to be included in the returned data.frame, 
#' specify them with \code{exclude}. 
#' @param symbols A vector of instrument names to include
#' @param exclude A vector of names of attributes that should not be included in the returned data.frame
#' @param attrs.of name of a FinancialInstrument instrument. Returned data.frame columns will be the attributes of instrument.
#' @return data.frame
#' @author Garrett See
#' @seealso \code{\link{buildHierarchy}}, \code{\link{instrument}}
#' @examples
#'
#' \dontrun{
#' currency('USD')
#' stock('GM','USD',exchange='NYSE')
#' stock('XOM','USD',description='Exxon Mobil')
#' instrument.table()
#' #Usually, currencies will not have as many attribute levels
#' #as other instruments, so you may want to exclude them from the table.
#' it <- instrument.table(exclude="USD|GM", attrs.of = "XOM") #columns created based on XOM instrument
#' #it <- instrument.table(exclude=c('USD','GM'), attrs.of = "XOM") #same thing
#' it <- instrument.table(exclude='tick_size|description|exchange')
#' }
#' @export
instrument.table <- function(symbols=NULL, exclude = NULL, attrs.of = NULL) {
#TODO check for numeric/character
    if (is.null(symbols)) symbols <- ls_instruments()
    if (is.null(attrs.of)) attrs.of <- symbols
   
    attr.names <- NULL
    for (symbol in attrs.of) {
        instr <- try(getInstrument(symbol,silent=TRUE),silent=TRUE)            
        if (!inherits(instr,'try-error') 
                && inherits(instr, 'instrument') ) {
            attr.names <- unique(c(attr.names, names(instr)))               
        }
    }
    if (length(exclude) > 1) 
        exclude <- paste(exclude, collapse='|')
    if (length(exclude == 1)) #i.e. if (!is.null(exclude))
        attr.names <- attr.names[-grep(exclude, attr.names, ignore.case=TRUE)]
    out <- buildHierarchy(symbols,attr.names)
  #FIXME: doesn't work if there is only 1 symbol.
    data.frame(out[,-1], row.names=as.character(out[,1]))
}


