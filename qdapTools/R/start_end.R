#' Get Location of Start/End Points
#' 
#' Get the locations of start/end places for the ones in a binary vector.
#' 
#' @param x A vector of 1 and 0 or \code{\link[base]{logical}}.
#' @return Returns a two column \code{\link[base]{data.frame}} of start and end locations for ones.
#' @references \url{http://stackoverflow.com/a/29184841/1000343}
#' @keywords start stop begin end
#' @export
#' @author Roland (\url{http://stackoverflow.com/users/1412059/roland}) and Tyler Rinker <tyler.rinker@@gmail.com>.
#' @examples
#' set.seed(10); (x <- sample(0:1, 50, TRUE, c(.35, .65)))
#' start_end(x)
#' (y <- sample(c(TRUE, FALSE), 50, TRUE, c(.35, .65)))
#' start_end(y)
start_end <- function(x){
    data.frame(
        start = which(diff(c(0L, x)) == 1L),
        end = which(diff(c(x, 0L)) == -1L)
    )
}
