#' Pacific Decadal Oscillation Index
#'
#' Monthly Pacific Decadal Oscillation (PDO) index
#' values from January 1900 to present.
#'
#' For more information see \url{https://github.com/poissonconsulting/rpdo}.
#'
#' @format A tbl data frame:
#' \describe{
#'   \item{Year}{The year as an integer.}
#'   \item{Month}{The month as an integer.}
#'   \item{PDO}{The Pacific Decadal Oscillation index as a numeric.}
#' }
#' @examples
#' library(rpdo)
#' library(ggplot2)
#'
#' data(pdo)
#' ggplot(data = subset(pdo, pdo$Month == 1), aes(x = Year, y = PDO)) +
#'  geom_line() + ylab("January PDO Index")
"pdo"
