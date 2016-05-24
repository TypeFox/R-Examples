#' 
#' Generate a data frame of pseudo-random numbers using the MLS Junk Generator algorithm
#' 
#' Based on user input seeds, this function generates a data frame of n 
#' pseudo-random numbers and names the column containing these as "RN" for "random
#' numbers."  This is achieved by calling junkgen.
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
#' # Generate a pseudo-random number data frame with 10 observations from user-specified seeds
#' 
#' w <- 1
#' x <- 2
#' y <- 3
#' z <- 4
#' 
#' mlsjunkgend(n = 10, w = w, x = x, y = y, z = z) # returns a data frame of 10 observations 
#' 
#' # Specifying different values for n and round
#' 
#' mlsjunkgend(n = 5, w = w, x = x, y = y, z = z, round = 2)
#' # returns a data frame identical to the above example but with only 5 observations 
#' # rounded to 2 decimal places
#'
#' # using the default value of n (1) is identical to assigning the rounded result of 
#' # junkgen to a data frame of 1 observation
#' 
#' round(junkgen(w = w, x = x, y = y, z = z), 5) # returns "[1] 0.95516"  
#' mlsjunkgend(w = w, x = x, y = y, z = z) 
#' # returns the following:
#' #        RN
#' # 1 0.95516
#' 
#' @return A numeric vector containing a single pseudo-random number

mlsjunkgend <- function(n = 1, w, x, y, z, round = 5) {
    if (is.numeric(n)) { mls <- data.frame()
                         for (i in 1:n) { ri <- junkgen(w, x, y, z)
                                          mls <- rbind(mls, round(ri, round))
                                          w <- x
                                          x <- y
                                          y <- z
                                          z <- ri
                                        } 
                         names(mls) <- "RN"
                         return(mls) 
                       }
    else { stop("Invalid input.  Please ensure n is numeric.") }
}