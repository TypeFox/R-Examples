#' Nonrandom year with sample
#' 
#' Nerdy way to wish someone a happy new year by using sample
#' 
#' @details Nerdy way to wish someone a happy new year, eg:\cr
#' Have a great \cr
#' \code{set.seed(1244); sample(0:9,4,T)}
#'
#' @return \code{\link{cat}}s command into the console that can be copypasted to anyone's R script.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, April 2014
#' @seealso \code{\link{nameSample}} to impress with "randomly" finding a name, 
#'          \code{\link{set.seed}}, \code{\link{sample}}, \code{\link{letters}}
#' @export
#' @examples
#' 
#' yearSample(2016)
#' # Have a nerdy
#' set.seed(12353); sample(0:9, 4, replace=TRUE)
#' 
#' @param year 4 digit numerical year.
#' 
yearSample <- function(
year
)
{
year2 <- as.numeric(substring(year, first=1:4, last=1:4))
year_is_false <- function(i)
  {
  set.seed(i)
  any(  sample(0:9, 4, replace=TRUE) != year2  )
  }
i <- 0
while( year_is_false(i) ) i <- i+1
output <- paste0("set.seed(", i, "); sample(0:9,4,T)\n")
message(output)
}

