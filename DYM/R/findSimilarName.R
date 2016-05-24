#' Looks for approximate matches to \code{x} (the first argument)
#' within \code{name} (the second) argument.
#' 
#' @param x
#'    A string giving the (misspelt) name to search for.
#' @param names
#'    A character vector of correct names to match against.
#' @param threshold
#'    The maximum distance between the misspell (\code{x}) and the correct answer (in \code{name}).
#' @param max_n
#'    An integer limiting the number of results.  Passed to \code{\link[utils]{head}}.
#' @param ignoreCase
#'    A logical value indicating whether differences in case should be ignored when matching.  Passed to \code{\link[utils]{adist}}.
#' @seealso \code{\link[utils]{adist}} calculates the distance between strings.
#' \code{\link[base]{agrep}} and \code{\link[stringdist]{stringdist-package}} 
#' provide alternate metrics for these distances.
#' @examples
#' \donttest{
#' x <- "logg"
#' names <- DYM:::getNames(x)
#' # Increasing threshold increases the number of hits, upto max_n = 10
#' lapply(
#'   stats::setNames(0:4, 0:4), 
#'   function(i) DYM:::findSimilarName(x, names, threshold = i)
#' )
#' 
#' # Use max_n = Inf to return all hits
#' DYM:::findSimilarName(x, names, threshold = 3, max_n = Inf)
#' 
#' # Negative max_n returns all hits except the last max_n
#' DYM:::findSimilarName(x, names, threshold = 3, max_n = -40)
#' 
#' # Set ignoreCase = TRUE to get more matches that differ by letter case
#' DYM:::findSimilarName(x, names, ignoreCase = TRUE)
#' }
#' @importFrom stats na.omit
#' @importFrom utils adist
#' @importFrom utils head
findSimilarName <- function (x, names, threshold = 2, max_n = 10, ignoreCase = FALSE) 
{
   if (is.na(x)) {
      character(0)
   }
   else {
      names <- na.omit(names)
      d <- as.integer(adist(x, names, ignore.case = ignoreCase))
      distance_data <- data.frame(
         names = names, 
         distance = d,
         stringsAsFactors = FALSE
      )[order(d), ]
      distance_data <- distance_data[distance_data$distance <= threshold, ]
      head(distance_data, max_n)$names
   }
}
