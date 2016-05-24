#' 
#' Generate a single pseudo-random number using the MLS Junk Generator algorithm
#' 
#' Based on user input seeds, this function generates a pseudo-random number.
#' This is called by the mlsjunkgen package's other functions to generate a 
#' pseudo-random number stream.
#' 
#' @param w the first seed required by the MLS Junk Generator algorithm
#' @param x the first seed required by the MLS Junk Generator algorithm
#' @param y the first seed required by the MLS Junk Generator algorithm
#' @param z the first seed required by the MLS Junk Generator algorithm
#' 
#' @export
#' 
#' @examples
#' # Generate a pseudo-random number with user-specified seeds
#' 
#' w <- 1
#' x <- 2
#' y <- 3
#' z <- 4
#' junkgen(w = w, x = x, y = y, z = z) # returns "[1] 0.9551644"
#' 
#' @return A numeric vector containing a single pseudo-random number

junkgen <- function (w, x, y, z) {
    if (is.numeric(w) & is.numeric(x) & is.numeric(y) & is.numeric(z)) {
        r <- 5.980217 * (w ^ 2) + 9.446377 * (x ^ 0.25) + 4.81379 * (y ^ 0.33) + 
         8.91197 * (z ^ 0.5)
        ri <- r - trunc(r)
    } else { stop("Invalid input.  Please ensure all seeds are numeric.") } 
    return(ri)
}