#' Compound Annual Growth Rate
#'
#' \itemize{
#'     \item \bold{geometric} FV = PV * (1 + geometric) ** years
#'     \item \bold{continuous} FV = PV * exp(continuous * years)
#'  }
#'
#' @note see \emph{r_continuous} and \emph{r_discrete}
#' @param PV the price at the beginning of the period
#' @param FV the price at the end of the period
#' @param fractional_years the length of the period in (fractional) years
#' @param type either "geometric" or "continuous"
#' @return the compounded rate of return, annualized
#'
#' @author George Fisher \email{GeorgeRFisher@gmail.com}
#'
#'
#' @examples
#' PV    <- 9000
#' FV    <- 13000
#' years <- 3
#' (geometric  <- CAGR(9000, 13000, years, type="geometric"))
#' (continuous <- CAGR(9000, 13000, years, type="continuous"))
#' 9000 * (1 + geometric) ** years
#' 9000 * exp(continuous * years)
#'
#' \dontrun{
#' error <- CAGR(9000, 13000, years, type="error")}
#'
#' @export
CAGR <- function(PV, FV, fractional_years, type="geometric") {
    if (type == "geometric")
        return (((FV / PV) ** (1 / fractional_years)) - 1)

    if (type == "continuous")
        return (log(FV / PV) / fractional_years)

    stop('CAGR takes type = ["geometric" | "continuous"]')
}
