#' 
#' Generate a vector of pseudo-random numbers using the MLS Junk Generator algorithm
#' 
#' Based on user input seeds, this function generates a vector of n 
#' pseudo-random numbers by calling junkgen.
#' 
#' @param n the number of pseudo-random numbers to generate; defaults to 1
#' @param w the first seed required by the MLS Junk Generator algorithm
#' @param x the first seed required by the MLS Junk Generator algorithm
#' @param y the first seed required by the MLS Junk Generator algorithm
#' @param z the first seed required by the MLS Junk Generator algorithm
#' @param round the number of decimal places to which to round the pseudo-random numbers; default = 5
#' 
#' @export
#' 
#' @examples
#' # Generate a pseudo-random number stream of length 5 with user-specified seeds
#' 
#' w <- 1
#' x <- 2
#' y <- 3
#' z <- 4
#'
#' # the following call returns "[1] 0.95516 0.66908 0.21235 0.34488 0.11995" 
#' mlsjunkgenv(n = 5, w = w, x = x, y = y, z = z) 
#' 
#' # Specifying different values for n and round
#' 
#' mlsjunkgenv(n = 3, w = w, x = x, y = y, z = z, round = 2) # returns "[1] 0.96 0.67 0.21"
#'
#' # using the default value of n (1) is identical to running junkgen and rounding 
#' # the result to 5 decimal places
#' 
#' round(junkgen(w = w, x = x, y = y, z = z),5) # returns "[1] 0.95516"  
#' mlsjunkgenv(w = w, x = x, y = y, z = z) # returns "[1] 0.95516"
#' 
#' @return A numeric vector containing a single pseudo-random number

mlsjunkgenv <- function(n = 1, w, x, y, z, round = 5) {
    if (is.numeric(n)) { mls <- numeric()
                         for (i in 1:n) { ri <- junkgen(w, x, y, z)
                                          mls <- c(mls, round(ri, round))
                                          w <- x
                                          x <- y
                                          y <- z
                                          z <- ri
                         } 
                         return(mls) 
                    }
    else { stop("Invalid input.  Please ensure n is numeric.") }
}