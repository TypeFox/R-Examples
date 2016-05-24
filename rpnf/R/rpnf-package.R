#' rpnf is a tool set to create and analyze Point & Figure Charts for given time series or data frame objects.
#' 
#' @name rpnf-package
#' @aliases rpnf-package
#' @aliases rpnf
#' @docType package
#' @title rpnf - The R Point & Figure Package
#' @author Sascha Herrmann \email{sascha.herrmann.consulting@@gmail.com}
#' @references Project Home Page \url{http://rpnf.r-forge.r-project.org}
#' @references Dorsey, Thomas J. Point and Figure Charting: The Essential Application for Forecasting and Tracking Market Prices. 3rd ed. Wiley Trading. Hoboken, N.J: John Wiley & Sons, 2007.
#' @references German version, which is the base for the package: Dorsey, Thomas. Sicher anlegen mit point & figure: klare Signale mit einfachen Methoden. Munich: FinanzBuch-Verl., 2000.
#' @keywords package
#' @seealso \code{\link{pnfprocessor}}
#' @seealso \code{\link{pnfplot}}
#' @seealso \code{\link{pnfplottxt}}
#' @examples
#' # Load rpnf library
#' library(rpnf) 
#' # Load free available sample data 
#' data(DOW) 
#' # Determine point and figure informations for a linear chart with boxsize of 1 point
#' pnfdata <- pnfprocessor(
#'   high=DOW$High,
#'   low=DOW$Low,
#'   date=DOW$Date,
#'   boxsize=1L,
#'   log=FALSE)  
#' # Show the object obtained
#' str(pnfdata)
#' # Show the data obtained
#' pnfdata
#' # Create a TXT based plot with X and O's
#' pnfplottxt(pnfdata,boxsize=1L,log=FALSE)
#' # Create a more graphical plot 
#' pnfplot(pnfdata)
#' \dontrun{
#' ### Second example: logarithmc example
#' # For most stocks and indices it is useful
#' # to do the analysis on a logarithmic scale.
#' # This can be done with pnfprocessor, too. 
#' # Ensure to make use of the getLogBoxsize() function 
#' # for an appropriate boxsize of a logarithmic chart.
#' # Determine point and figure informations for a logarithmic chart with boxsize 2\% 
#' symbol.pnf <- pnfprocessor(
#'   high=DOW$High,
#'   low=DOW$Low,
#'   date=DOW$Date,
#'   boxsize=getLogBoxsize(2),
#'   log=TRUE)  
#' 
#' # View the result
#' tail(symbol.pnf)
#' #View(symbol.pnf)
#' 
#' # or plot it as a modern chart
#' pnfplot(symbol.pnf,main="P&F Plot DOW (log)")
#' # Or in the old traditional TXT style
#' pnfplottxt(symbol.pnf,boxsize=getLogBoxsize(2),log=TRUE,main="P&F Plot DOW (log)")
#' 
#' ### Additional examples
#' # Examples for additional uses cases like
#' # - relative strength vs index
#' # - bullish percent of an index
#' # - and many others 
#' # can be found in your local package library directory.
#' # Search for rpnf-example1.R, rpnf-example2.R and so on.
#' }
NULL

#' This is some free available quote data for the DOW Chemical Company.
#' 
#'  End of day open, high, low, close and volume, dividends and splits, 
#'  and split/dividend adjusted open, high, low close and volume for Dow Chemical Company (The) (DOW).  
#'  Data are freely available at https://www.quandl.com/data/WIKI/DOW,
#'  and may be copy, distribute, disseminate or include the data in other products for commercial and/or noncommercial purposes. 
#'  This data is part of Quandl's Wiki initiative to get financial data permanently into the public domain. 
#'  Quandl relies on users like you to flag errors and provide data where data is wrong or missing. 
#'  Get involved: connect@@quandl.com 
#'
#' @name DOW
#' @docType data
#' @author Sascha Herrmann \email{sascha.herrmann.consulting@@gmail.com}
#' @references \url{https://www.quandl.com/data/WIKI/DOW}
#' @keywords data
NULL
