#' Simple Random Sample of Golf Courses
#' 
#' Simple Random Sample (SRS) of 120 golf courses taken from 
#'  the population of the Website www.golfcourse.com
#' @name golfsrs
#' @docType data
#' @format Data frame with the following 16 variables: 
#' \describe{
#'   \item{RN}{random number used to select golf course for sample}
#'   \item{state}{state name}
#'   \item{holes}{number of holes}
#'   \item{type}{type of course; factor with levels \code{priv} (private),
#'     \code{semi} (semi-private), \code{pub} (public), \code{mili}
#'     (military) and \code{res} (resort)}
#'   \item{yearblt}{year the course was built}
#'   \item{wkday18}{greens fee for 18 holes during week}
#'   \item{wkday9}{greens fee for 9 holes during week}
#'   \item{wkend18}{greens fee for 18 holes on weekend}
#'   \item{wkend9}{greens fee for 9 holes on weekend}
#'   \item{backtee}{back-tee yardage}
#'   \item{rating}{course rating}
#'   \item{par}{par for course}
#'   \item{cart18}{golf cart rental fee for 18 holes}
#'   \item{cart9}{golf cart rental fee for 9 holes}
#'   \item{caddy}{Are caddies available? factor with levels \code{yes}
#'     and \code{no}}
#'   \item{pro}{Is a golf pro available? factor with levels \code{yes}
#'     and \code{no}}
#' }
#' @source \url{http://www.golfcourse.com}
#' @references Lohr (1999). Sampling: Design and Analysis, Duxbury, p. TODO and
#'   TODO.
#' @export
NULL
