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
# $Id: buildHierarchy.R 963 2012-03-04 21:01:36Z gsee $
#
###############################################################################


#' Construct a hierarchy of instruments useful for aggregation
#'
#' Construct a hierarchy of instruments useful for aggregation
#' 
#' @param primary_ids A character vector of \code{instrument} primary_ids to be 
#' included in the hierarchy list
#' @param ... character names of instrument attributes in top-down order. 
#' @author Peter Carl, Alexis Petit, Garrett See
#' @return Constructs a data.frame that contains the list of assets in the first 
#' column and the category or factor for grouping at each level in the following 
#' columns
#' @seealso \code{\link{instrument.table}}
# TODO add a link to PortfolioAnalytics attribution functions, when they exist
#' @export
#' @examples
#' \dontrun{
#' # rm_instruments(keep.currencies=FALSE)
#' ## Define some stocks
#' update_instruments.TTR(c("XOM", "IBM", "CVX", "WMT", "GE"), exchange="NYSE")
#' 
#' buildHierarchy(ls_instruments(), "type")
#' buildHierarchy(ls_stocks(), c("Name", "Sector"))
#' buildHierarchy(ls_stocks(), "Industry", "MarketCap")
#' }
buildHierarchy <- function(primary_ids, ...) {
    levels <- unlist(list(...))
    if (!is.null(levels)) stopifnot(is.character(levels))
    out <- data.frame(primary_ids, stringsAsFactors=FALSE)
    ilist <- lapply(primary_ids, getInstrument)
    for (level in levels) {
        tmp_level <- as.character(lapply(1:length(primary_ids), 
                                         function(x) ilist[[x]][[level]]))
        out <- cbind(out, tmp_level, stringsAsFactors=FALSE)
    }
    colnames(out) <- c("primary_id", levels)
    return(out)
}
