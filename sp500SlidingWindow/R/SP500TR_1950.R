#' Daily S&P 500 Total Return data from Jan 3, 1950 to present
#'
#' @description
#' Yahoo Finance returns TR data (^SP500TR) from 1988, non-TR (^GSPC) data
#' from 1950. I spent a lot of time working with Bob Schiller's data
#' (\url{http://www.econ.yale.edu/~shiller/data.htm}) which contains dividends
#' back to 1871. What I found was that adjusting for dividends was difficult
#' and I could never improve upon appending the TR data to the non-TR data
#' at the 1987-1988 year break. Neither can anyone else I can find on the Internet.
#' I keep hoping to find better data than this but so far I have been stymied.
#'
#' @return A data.frame with the daily data
#'
#' @references
#' Yahoo Finance, Bob Schiller
#'
#' @details The columns of the data.frame returned
#' \itemize{
#'     \item \bold{Date}
#'     \item \bold{Open}
#'     \item \bold{High}
#'     \item \bold{Low}
#'     \item \bold{Close}
#'     \item \bold{Volume}
#'     \item \bold{Adj.Close}
#'     \item \bold{Year}
#'     \item \bold{Month}
#'     }
#'
#' @importFrom utils read.table
#'
#' @author George Fisher
#'
#' @examples
#' sp500_idx <- SP500TR_1950()
#' head(sp500_idx)
#' tail(sp500_idx)
#'
#' @export
SP500TR_1950 <- function() {

    Date <- NULL

    # ================== read the 1988 - Present TR data ==================

    URL_tr <- "http://ichart.finance.yahoo.com/table.csv?s=%5ESP500TR"

    SP500TR <- read.table(URL_tr,
                          header = TRUE, sep=",", stringsAsFactors = FALSE)
    SP500TR$Date <- as.Date(SP500TR$Date, "%Y-%m-%d")
    SP500TR     <- dplyr::arrange(SP500TR, Date)

    SP500TR$Year  <- lubridate::year(SP500TR$Date)
    SP500TR$Month <- lubridate::month(SP500TR$Date)

    # ================== read the 1950 - Present non-TR data ==================

    URL_sp <- "http://ichart.finance.yahoo.com/table.csv?s=%5EGSPC"

    GSPC <- read.table(URL_sp,
                       header = TRUE, sep=",", stringsAsFactors = FALSE)
    GSPC$Date <- as.Date(GSPC$Date, "%Y-%m-%d")
    GSPC     <- dplyr::arrange(GSPC, Date)

    GSPC$Year  <- lubridate::year(GSPC$Date)
    GSPC$Month <- lubridate::month(GSPC$Date)

    # ================== combine ==================
    # ======= non-TR 1950-1987 TR 1988-Present ====

    tr <- rbind(GSPC[GSPC$Year<min(SP500TR$Year),], SP500TR)

    return(tr)
}
