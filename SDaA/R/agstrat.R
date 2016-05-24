#' Data from a stratified random sample of size 300 from 
#' the U.S. 1992 Census of Agriculture.
#' 
#' @name agstrat
#' @docType data
#' @format Data frame with the following 17 variables: 
#' \describe{
#'   \item{county}{county name}
#'   \item{state}{state abbreviation}
#'   \item{acres92}{number of acres devoted to farms, 1992}
#'   \item{acres87}{number of acres devoted to farms, 1987}
#'   \item{acres82}{number of acres devoted to farms, 1982}
#'   \item{farms92}{number of farms, 1992}
#'   \item{farms87}{number of farms, 1987}
#'   \item{farms82}{number of farms, 1982}
#'   \item{largef92}{number of farms with 1000 acres or more, 1992}
#'   \item{largef87}{number of farms with 1000 acres or more, 1987}
#'   \item{largef82}{number of farms with 1000 acres or more, 1982}
#'   \item{smallf92}{number of farms with 9 acres or fewer, 1992}
#'   \item{smallf87}{number of farms with 9 acres or fewer, 1987}
#'   \item{smallf82}{number of farms with 9 acres or fewer, 1982}
#'   \item{region}{factor with levels \code{S} (south), \code{W} (west), 
#'     \code{NC} (north central), \code{NE} (northeast)}
#'   \item{rn}{random numbers used to select sample in each stratum}
#'   \item{weight}{sampling weighs for each county in sample}
#' }
#' @source U.S. 1992 Census of Agriculture
#' @references Lohr (1999). Sampling: Design and Analysis, Duxbury, p. TODO and
#'   437.
#' @export
NULL
