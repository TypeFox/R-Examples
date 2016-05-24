#' Shift Vector Left/Right 
#' 
#' Shift a vector left or right n spaces.
#' 
#' @param x A vector.
#' @param n The number of moves left or right to shift.
#' @param direction A direction to shift; must be either "left" or "right".
#' Use explicit directional shift functions \code{shift_right} and 
#' \code{shift_left} for better performance.
#' @return Returns a shifted vector.
#' @rdname shift
#' @keywords shift
#' @export
#' @rdname shift
#' @examples
#' lapply(0:9, function(i) shift(1:10, i))
#' lapply(0:9, function(i) shift(1:10, i, "left"))
#' 
#' ## Explicit, faster shifting
#' lapply(0:9, function(i) shift_right(1:10, i))
#' lapply(0:9, function(i) shift_left(1:10, i))
#' lapply(0:25, function(i) shift_left(LETTERS, i))
shift <- function(x, n, direction = "right") {

    if (!direction %in% c("left", "right")) {
        stop("`direction` must be either \"left\" or \"right\"")
    }
    match.fun(paste0("shift_", direction))(x=x, n=n)
}

#' @export
#' @rdname shift
shift_right <- function(x, n){
    if (n == 0) return(x)
    c(x[(n+1):length(x)], x[1:n])
}

#' @export
#' @rdname shift
shift_left <- function(x, n){
    if (n == 0) return(x)
    c(x[(length(x)-(n-1)):length(x)], x[1:(length(x)-n)])
}
