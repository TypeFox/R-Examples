#' Generate all point and figure informations for a given time series.
#' 
#' Please ensure that high, low and date are all ordered according to the Date column.
#' 
#' @param high a vector containing the high quotes
#' @param low a (optional) vector containing the low quotes
#' @param date a vector of dates the quotes belong
#' @param reversal number of boxes needed to make a reversal 
#' @param boxsize the boxsize to be used
#' @param log should we do the calculations on a logarithmic scale
#' @param style the style the pnfprocessor is working with. Can be \{xo,rs,bp\}.
#' @return returns a data table with all point and figure information in it
#' @export
#' @seealso \code{\link{pnfplot}}
#' @seealso \code{\link{pnfplottxt}}
#' @references \url{http://rpnf.r-forge.r-project.org}
#' @examples
#' library(rpnf) # Load rpnf library
#' data(DOW) # (Offline) Load free available sample data from https://www.quandl.com/data/WIKI/DOW
#' pnfdata <- pnfprocessor(
#'   high=DOW$High,
#'   low=DOW$Low,
#'   date=DOW$Date,
#'   boxsize=1L,
#'   log=FALSE)  
#' pnfdata
pnfprocessor <- function(
  high,
  low=high,
  date,
  reversal=3L, 
  boxsize=1L, 
  log=FALSE,
  style="xo") {
  # check for proper style selection
  if (!style %in% c("xo","rs","bp")) {
    stop("Select a proper chart style: 'xo', 'rs' or'bp'!")
  }
  # first execute basic time series processing
  result <- xo.processor(high=high,low=low,date=date,reversal=reversal,boxsize=boxsize,log=log)
  # now check for style to enhance the result
  if (style == "xo") {
    result <- xo.signalprocessor(result,reversal) # 1.2 sec
    result <- xo.trendline.processor(result) # 0.25 sec
    result <- xo.priceobjective.processor(result,reversal,boxsize,log) # 2.2 sec
  } else if (style == "bp") {
    result <- bp.signalprocessor(result)
  } else if (style == "rs") {
    result <- rs.signal.processor(result)
  } else {
    warning("Unknown style detected! No chart enhancements were made!")
  }
  result
}

