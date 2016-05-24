#' Lowest nonzero values
#' 
#' Sometimes want to plot p-values (e.g., volcano plot or MA-plot), but if a
#' statistical test returns a zero p-value, this causes problems with
#' visualization on the log scale. This function returns a vector where  the
#' zero values are equal to the smallest nonzero value in the vector.
#' 
#' @author Stephen Turner
#' @keywords keywords
#'   
#' @param x A vector of p-values between 0 and 1.
#'   
#' @return A vector of p-values where zero values are exchanged for the lowest non-zero p-value in the original vector.
#' 
#' @importFrom stats na.omit
#'   
#' @examples
#' lowestnonzero(c(.042, .02, 0, .001, 0, .89))
#' 
#' @export
lowestnonzero <- function(x) {
    if(any(na.omit(x)<0) | any(na.omit(x)>1)) stop("P-values should be between 0 and 1.")
    xo <- x[order(x)]
    lnz <- xo[xo>0][1]
    x[x==0] <- lnz
    x
}
