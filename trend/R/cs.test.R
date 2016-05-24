cs.test <- function(x) {
##    Copyright (C) 2015, 2016  Thorsten Pohlert
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##    This function computes the Cox Stuart trend test.
##
    dname <- deparse(substitute(x))
    stopifnot(is.numeric(x))
    x <- na.omit(x)
    n <- length(x)
    if(n < 3) {
        stop("sample size must be greater than 2")
    }
    l <- ceiling(n/3)
    u <- x[(n-l+1):n] - x[1:l]
    u <- u[u != 0]  # remove 0 values
    if (length(u) == 0){
        stop("entire sample contains identical values") 
    }
    sgn <- sign(u)
    S <- table(sgn)
    zstat <- function(S, n) {
        if (n > 30) {
            out <- abs(S - n / 6) / sqrt(n / 12)
            return(out)
        } else {
            out <- (abs(S - n / 6) - 0.5) / sqrt(n / 12)
            return(out)
        }
    }
    S <- max(S)
    z <- zstat(S, n)
    if (z >= 0) {
        p.value <- 2 * pnorm(z, lower.tail=FALSE)
    } else {
        p.value <- 2 * pnorm(z, lower.tail=TRUE)
    }
   
    names(z) <- "Cox and Stuart z-value"
    out <- list(statistic = z, p.value = p.value,
                alternative = "monotonic trend",
                data.name = dname,
                method = "Cox and Stuart trend test")
    class(out) <- "htest"
    out
}
