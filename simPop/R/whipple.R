#' Whipple index (original and modified)
#' 
#' The function calculates the original and modified Whipple index to evaluate
#' age heaping.
#' 
#' The original Whipple's index is obtained by summing the number of persons in
#' the age range between 23 and 62, and calculating the ratio of reported ages
#' ending in 0 or 5 to one-fifth of the total sample. A linear decrease in the
#' number of persons of each age within the age range is assumed. Therefore,
#' low ages (0-22 years) and high ages (63 years and above) are excluded from
#' analysis since this assumption is not plausible.
#' 
#' When the digits 0 and 5 are not reported in the data, the original Whipple
#' index varies between 0 and 100, 100 if no preference for 0 or 5 is within
#' the data. When only the digits 0 and 5 are reported in the data it reaches a
#' to a maximum of 500.
#' 
#' For the modified Whipple index, age heaping is calculated for all ten digits
#' (0-9). For each digit, the degree of preference or avoidance can be
#' determined for certain ranges of ages, and the modified Whipple index then
#' is given by the absolute sum of these (indices - 1).
#' 
#' @name whipple
#' @param x numeric vector holding the age of persons
#' @param method \dQuote{original} or \dQuote{modified} Whipple index.
#' @return The original or modified Whipple index.
#' @author Matthias Templ
#' @seealso \code{\link{sprague}}
#' @references Henry S. Shryock and Jacob S. Siegel, Methods and Materials of
#' Demography (New York: Academic Press, 1976)
#' @keywords arith
#' @export
#' @examples
#' 
#' age <- sample(1:100, 5000, replace=TRUE)
#' whipple(age)
#' 
whipple <- function(x, method="standard"){
  x <- x[x >= 23 & x <= 62]
  if(method == "standard"){
    xm <- x %% 5
    whipple <- (length(xm[xm==0])/length(x))*500
  }
  if(method == "modified"){
    tab <- table(x)
    sp <- function(p) seq(p,p+30,10)
    W <- numeric(9)
    W[1] <- 5*sum(sp(31)) / sum(5*sp(29)) 
    W[2] <- 5*sum(sp(32)) / sum(5*sp(30)) 
    W[3] <- 5*sum(sp(23)) / sum(5*sp(21))
    W[4] <- 5*sum(sp(24)) / sum(5*sp(22))
    W[6] <- 5*sum(sp(26)) / sum(5*sp(24))
    W[7] <- 5*sum(sp(27)) / sum(5*sp(25))
    W[8] <- 5*sum(sp(28)) / sum(5*sp(26))
    W[9] <- 5*sum(sp(29)) / sum(5*sp(27))
    whipple <- sum(abs(W-1), na.rm=TRUE)
  }
  return(whipple)
}

