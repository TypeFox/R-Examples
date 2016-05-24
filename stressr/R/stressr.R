#' @title Get financial stress index data.
#' @description Downloads Cleveland FRB financial stress index data.
#' @details  Transforms the HTML into a data frame, transforms the character date into \code{Date} objects, and then an \code{xts} object.  
#' @param verbose whether to print progress messages, default FALSE
#' @return List of class type \code{cfsi} containing \code{xts} time history object \code{df}, plot colors array \code{colors}, default plot main title \code{main}, and default plot y-axis label \code{ylab}.
#' @export
#' @importFrom XML readHTMLTable
#' @import xts
#' @references \url{http://www.clevelandfed.org/research/data/financial_stress_index/index.cfm}
#' @seealso getStressComponents
#' @examples
#' \dontrun{ 
#' getStressIndex()
#' }
#' 
getStressIndex <- function(verbose=FALSE) {
  
  ciurl <- "http://www.clevelandfed.org/research/data/financial_stress_index/cfsi-dataset.cfm"
  
  if ( verbose )
    message("Reading CFSI data...")
  x <- readHTMLTable(ciurl,
                     header=TRUE,
                     skip.rows=1,
                     colClasses=c("character",rep("numeric",1)),
                     stringsAsFactors=FALSE)
  
  cfsi <- x[[1]]
  colnames(cfsi) <- c(
    "DATE",
    "CLEVELAND")
  
  cfsi$DATE <- as.Date(cfsi$DATE,format="%m/%d/%Y")
  cfsi <- xts(cfsi[,-1],order.by=cfsi[,1])
  colnames(cfsi) <- "CLEVELAND"
  
  if (verbose)
    message("Completed CFSI data.")
  
  rv <- list(df=cfsi,
             colors="#376e90",
             main="Financial Stress Index",
             ylab="(Z-Score)")
  class(rv) <- "cfsi"
  rv
}

#' @title Get financial stress index component data.
#' @description Downloads Cleveland FRB financial stress index data.
#' @details  Transforms the HTML into a data frame, transforms the character date into \code{Date} objects, and then an \code{xts} object.  
#' @param verbose whether to print progress messages, default FALSE
#' @return List of class type \code{stress} containing \code{xts} time history object \code{df}, plot colors array \code{colors}, default plot main title \code{main}, and default plot y-axis label \code{ylab}.
#' @export
#' @importFrom XML readHTMLTable
#' @import xts
#' @references \url{http://www.clevelandfed.org/research/data/financial_stress_index/index.cfm}
#' @seealso getStressIndex
#' @note Meant for internal use by the other, more specific, query functions.
#' @examples
#' \dontrun{ 
#' getStressComponents()
#' }
#' 
getStressComponents <- function(verbose=FALSE) {
  
  cfurl <- "http://www.clevelandfed.org/research/data/financial_stress_index/component-dataset.cfm"
  
  if (verbose)
    message("Reading component data...")
  x <- readHTMLTable(cfurl,
                     header=TRUE,
                     skip.rows=1,
                     colClasses=c("character",rep("numeric",22)),
                     stringsAsFactors=FALSE)
  components <- x[[1]]
  colnames(components) <- c(
    "DATE",
    "ASSET-BACKED SECURITY SPREAD",
    "BANK BOND SPREAD",
    "COMMERCIAL MBS SPREAD",
    "COMMERCIAL PAPER - T-BILL SPREAD",
    "COMMERCIAL REAL ESTATE SPREAD",
    "CORPORATE BOND SPREAD",
    "COVERED INTEREST SPREAD",
    "CREDIT MARKETS",
    "EQUITY MARKETS",
    "FINANCIAL BETA",
    "FOREIGN EXCHANGE MARKETS",
    "INTERBANK COST OF BORROWING",
    "INTERBANK LIQUIDITY SPREAD",
    "INTERBANK MARKETS",
    "LIQUIDITY SPREAD",
    "REAL ESTATE MARKET",
    "RESIDENTIAL MBS SPREAD",
    "RESIDENTIAL REAL ESTATE SPREAD",
    "SECURITIZATION MARKET",
    "STOCK MARKET CRASHES",
    "TREASURY YIELD CURVE SPREAD",
    "WEIGHTED DOLLAR CRASHES"
  )
  
  fedColors <- c(
    "#c1ad97", # "ASSET-BACKED SECURITY SPREAD",
    "#e4aa71", # "BANK BOND SPREAD",
    "#a8927a", # "COMMERCIAL MBS SPREAD",
    "#537f9a", # "COMMERCIAL PAPER - T-BILL SPREAD",
    "#597b6a", # "COMMERCIAL REAL ESTATE SPREAD",
    "#d0e7f7",# "CORPORATE BOND SPREAD",
    "#40657c", # "COVERED INTEREST SPREAD",
    "#5e8aa5", # "CREDIT MARKETS",
    "#dad38d", # "EQUITY MARKETS",
    "#bc8146", # "FINANCIAL BETA",
    "#999999", # "FOREIGN EXCHANGE MARKETS",
    "#d0945a", # "INTERBANK COST OF BORROWING",
    "#fac38e", #"INTERBANK LIQUIDITY SPREAD",
    "#a0c81d", # "INTERBANK MARKETS",
    "#324754", # "LIQUIDITY SPREAD",
    "#2ba4a5", # "REAL ESTATE MARKET",
    "#97846e", # "RESIDENTIAL MBS SPREAD",
    "#6e9180", # "RESIDENTIAL REAL ESTATE SPREAD",
    "#df6968", #" SECURITIZATION MARKET",
    "#dbd48e", # "STOCK MARKET CRASHES",
    "#90b7cf", # "TREASURY YIELD CURVE SPREAD",
    "#999999" # WEIGHTED DOLLAR CRASHES"
  )

  
  # convert date posted to date class
  components$DATE <- as.Date(components$DATE,format="%m/%d/%Y")
  components <- xts(components[,-1],order.by=components[,1])
  if ( verbose )
    message("Completed component data.")
  
  rv <- list(df=components,
             colors=fedColors,
             columns=1:ncol(components),
             main="Financial Stress Index",
             ylab="")
  class(rv) <- "stress"
  rv
}

#' @title Get stress components summary
#' @description Downloads FRB financial stress index component data.
#' @details Downloads the Cleveland FRB data products for financial stress index components daily time series.  Component values include
#' \itemize{
#' \item foreign exchange markets
#' \item credit markets
#' \item interbank markets
#' \item equity markets
#' \item real estate market
#' \item securitization market
#' }
#' @param s the list of class \code{stress} from previous queries, or NULL to perform new query
#' @return A list of class \code{stress}
#' @export
#' @seealso getStressData getEquityMarkets getFundingMarkets getCreditMarkets getForeignExchangeMarkets getRealEstateMarkets getSecuritizationMarkets
#' @examples
#' \dontrun{
#' getEquityMarkets()
#' }
getComponentSummary <- function(s=NULL) {
  stopifnot(class(s)=="stress" || is.null(s))
  
    fedColors <- c(
      "#999999", # "FOREIGN EXCHANGE MARKETS"
      "#5e8aa5", # "CREDIT MARKETS"
      "#a0c81d", # "INTERBANK MARKETS"
      "#dad38d", # "EQUITY MARKETS"
      "#2ba4a5", # "REAL ESTATE MARKET"
      "#df6968"  #" SECURITIZATION MARKET"
    )
  
  if ( is.null(s) )
    s <- getStressComponents()
  
  s$columns = c(11,8,14,9,16,19)
  s$colors = fedColors
  s$main="Components of the CFSI - Summary"
  s$ylab=""
  s
}

#' @title Get equity markets stress components
#' @description Downloads FRB financial stress index component data.
#' @details Downloads the Cleveland FRB data products for financial stress index components daily time series.  Component values include
#' \itemize{
#' \item stock market crashes
#' }
#' @param s the list of class \code{stress} from previous queries, or NULL to perform new query
#' @return A list of class \code{stress}
#' @export
#' @seealso getStressData getEquityMarkets getFundingMarkets getCreditMarkets getForeignExchangeMarkets getRealEstateMarkets getSecuritizationMarkets
#' @examples
#' \dontrun{
#' getEquityMarkets()
#' }
getEquityMarkets <- function(s=NULL) {
  stopifnot(class(s)=="stress" || is.null(s))
  
  fedColors = c(
    "#dbd48e" # stock market crashes
  )
  
  if ( is.null(s) )
    s <- getStressComponents()
  
  s$columns = 9
  s$colors = fedColors
  s$main="Equity Markets"
  s$ylab=""
  s
}

#' @title Get funding markets stress components
#' @description Downloads FRB financial stress index component data.
#' @details Downloads the Cleveland FRB data products for financial stress index components daily time series.  Component values include
#' \itemize{
#' \item financial beta
#' \item interbank cost of borrowing
#' \item bank bond spread
#' \item interbank liquidity spread
#' }
#' @param s the list of class \code{stress} from previous queries, or NULL to perform new query
#' @return A list of class \code{stress}
#' @export
#' @seealso getStressData getEquityMarkets getCreditMarkets getForeignExchangeMarkets getRealEstateMarkets getSecuritizationMarkets
#' @examples
#' \dontrun{
#' getFundingMarkets()
#' }
getFundingMarkets <- function(s=NULL) {
  stopifnot(class(s)=="stress" || is.null(s))
  
  fedColors = c(
    "#bc8146", # financial beta,
    "#d0945a", # interbank cost of borrowing
    "#e4aa71", # bank bond spread
    "#fac38e" # interbank liquidity spread
  )
  
  if ( is.null(s) )
    s <- getStressComponents()
  
  s$columns = c(10,12,2,13)
  s$colors = fedColors
  s$main="Funding Markets"
  s$ylab=""
  s
}

#' @title Get credit markets stress components
#' @description Downloads FRB financial stress index component data.
#' @details Downloads the Cleveland FRB data products for financial stress index components daily time series.  Component values include
#' \itemize{
#' \item liquidity spread
#' \item covered interest spread
#' \item commercial paper - t-bill spread
#' \item treasury yield curve spread
#' \item coporate bond spread
#' }
#' @param s the list of class \code{stress} from previous queries, or NULL to perform new query
#' @return A list of class \code{stress}
#' @export
#' @seealso getStressData getEquityMarkets getFundingMarkets getForeignExchangeMarkets getRealEstateMarkets getSecuritizationMarkets
#' @examples
#' \dontrun{
#' getCreditMarkets()
#' }
getCreditMarkets <- function(s=NULL) {
  stopifnot(class(s)=="stress" || is.null(s))
  
  fedColors = c(
    "#324754", # liquidity spread
    "#40657c", # covered interest spread
    "#537f9a", # commercial paper - t-bill spread
    "#90b7cf", # treasury yield curve spread
    "#d0e7f7"  # corporate bond spread
  )
  
  if ( is.null(s) )
    s <- getStressComponents()
  
  s$columns = c(15,7,4,21,6)
  s$colors = fedColors
  s$main="Credit Markets"
  s$ylab=""
  s
}

#' @title Get foreign exchange markets stress components
#' @description Downloads FRB financial stress index component data.
#' @details Downloads the Cleveland FRB data products for financial stress index components daily time series.  Component values include
#' \itemize{
#' \item weighted dollar crashes
#' }
#' @param s the list of class \code{stress} from previous queries, or NULL to perform new query
#' @return A list of class \code{stress}
#' @export
#' @seealso getStressData getEquityMarkets getCreditMarkets getFundingMarkets  getRealEstateMarkets getSecuritizationMarkets
#' @examples
#' \dontrun{
#' getForeignExchangeMarkets()
#' }
getForeignExchangeMarkets <- function(s=NULL) {
  stopifnot(class(s)=="stress" || is.null(s))
  
  fedColors = c(
    "#999999" # weighted dollar crashes
  )
  
  if ( is.null(s) )
    s <- getStressComponents()
  
  s$columns = 22
  s$colors = fedColors
  s$main="Foreign Exchange Markets"
  s$ylab=""
  s
}


#' @title Get foreign exchange markets stress components
#' @description Downloads FRB financial stress index component data.
#' @details Downloads the Cleveland FRB data products for financial stress index components daily time series.  Component values include
#' \itemize{
#' \item commecial real estate spread
#' \item residential real estate spread
#' }
#' @param s the list of class \code{stress} from previous queries, or NULL to perform new query
#' @return A list of class \code{stress}
#' @export
#' @seealso getStressData getEquityMarkets getCreditMarkets getFundingMarkets  getForeignExchangeMarkets getSecuritizationMarkets
#' @examples
#' \dontrun{
#' getRealEstateMarkets()
#' }
getRealEstateMarkets <- function(s=NULL) {
  stopifnot(class(s)=="stress" || is.null(s))
  
  fedColors = c(
    "#597b6a", # commercial real estate spread
    "#6e9180" # residential real estate spread
  )
  
  if ( is.null(s) )
    s <- getStressComponents()
  
  s$columns = c(5,18)
  s$colors = fedColors
  s$main="Real Estate Markets"
  s$ylab=""
  s
}

#' @title Get securitization markets stress components
#' @description Downloads FRB financial stress index component data.
#' @details Downloads the Cleveland FRB data products for financial stress index components daily time series.  Component values include
#' \itemize{
#' \item residential MBS spread
#' \item commercial MBS spread
#' \item asset-backed security spread
#' }
#' @param s the list of class \code{stress} from previous queries, or NULL to perform new query
#' @return A list of class \code{stress}
#' @export
#' @seealso getStressData getEquityMarkets getCreditMarkets getFundingMarkets  getForeignExchangeMarkets getRealEstateMarkets
#' @examples
#' \dontrun{
#' getSecuritizationMarkets()
#' }
getSecuritizationMarkets <- function(s=NULL) {
  stopifnot(class(s)=="stress" || is.null(s))
  
  fedColors = c(
    "#97846e", # residential MBS spread
    "#a8927a", # commercial MBS spread
    "#c1ad97"  # assets-backed security spread
  )
  
  if ( is.null(s) )
    s <- getStressComponents()
  
  s$columns = c(17,3,1)
  s$colors = fedColors
  s$main="Securitization Markets"
  s$ylab=""
  s
}

