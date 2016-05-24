sens.slope <-
function(x, level=0.95){
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
##    This function computes Sens's Slope.
##
    na.fail(x)
    res.mk <- mk.test(x)
    varS <- res.mk$Varianz
    n <- length(x)
 #   d <- NULL
    retval <- list(b.sen = NULL, b.sen.up = NULL, b.sen.lo = NULL,
                   intercept = NULL, nobs = n, method="SSLP", D=NULL,
                   conf.int = level * 100, varS = NULL)
    k <- 0
    d <- rep(NA, n * (n-1)/2)
    for (i in 1:(n-1)) {
        for (j in (i+1):n){
            k <- k + 1
            d[k] <- (x[j] - x[i]) / ( j - i)
        }
    }
    b.sen <- median(d, na.rm=TRUE)
    C <- qnorm(1 - (1 - level)/2) * sqrt(varS)
    rank.up <- round((k + C) / 2 + 1)
    rank.lo <- round((k - C) / 2)
    rank.d <- sort(d)
    retval$b.sen.lo <- rank.d[rank.lo]
    retval$b.sen.up <- rank.d[rank.up]
    intercepts <- x - b.sen * (1:n)
    retval$intercept <- median(intercepts)
    retval$b.sen <- b.sen
    retval$D <- d
    retval$varS <- varS
    class(retval) <- "trend.test"
    return(retval)
}

