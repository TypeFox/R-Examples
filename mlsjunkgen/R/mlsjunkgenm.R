#' 
#' Generate a matrix of pseudo-random numbers using the MLS Junk Generator algorithm
#' 
#' Based on user input seeds, this function generates a vector of n 
#' pseudo-random numbers by calling mlsjunkgenv which in turn calls junkgen. 
#' 
#' @param nrow the number of rows for the matrix; defaults to 1
#' @param ncol the number of columns for the matrix; defaults to 1
#' @param w the first seed required by the MLS Junk Generator algorithm
#' @param x the first seed required by the MLS Junk Generator algorithm
#' @param y the first seed required by the MLS Junk Generator algorithm
#' @param z the first seed required by the MLS Junk Generator algorithm
#' @param round the number of decimal places to which to round the pseudo-random numbers; default = 5
#' 
#' @export
#' 
#' @examples
#' # Generate a 4x4 matrix of pseudo-random numbers with user-specified seeds
#' 
#' w <- 1
#' x <- 2
#' y <- 3
#' z <- 4
#' 
#' mlsjunkgenm(nrow = 4, ncol = 4, w = w, x = x, y = y, z = z) # returns a 4x4 matrix
#' 
#' # the sixteen values in the above matrix are equivalent to the following call 
#' # to mlsjunkgenv
#'  
#' mlsjunkgenv(n = 16, w = w, x = x, y = y, z = z)
#' 
#' # matrices need not be square
#' # this returns a 3x2 matrix of pseudo-random numbers with 2 decimal places
#' mlsjunkgenm(nrow = 3, ncol = 2, w = w, x = x, y = y, z = z, round = 2) 
#'
#' # using the default value of n (1) generates a 1x1 matrix the value of which 
#' # is identical to running junkgen and rounding the result to 5 decimal places
#' 
#' round(junkgen(w = w, x = x, y = y, z = z), 5) # returns "[1] 0.95516"  
#' mlsjunkgenv(w = w, x = x, y = y, z = z) # returns a 1x1 matrix with single element = "0.95516"
#' 
#' @return A numeric vector containing a single pseudo-random number

mlsjunkgenm <- function(nrow = 1, ncol = 1, w, x, y, z, round = 5) {
    if (is.numeric(nrow) & is.numeric(ncol)) { 
        matrix(mlsjunkgenv(nrow * ncol, w, x, y, z, round), nrow = nrow, 
               ncol = ncol)
    } else { stop("Invalid input.  Please ensure nrow and ncol are numeric.") }
}