sea.sens.slope <-
    function(x){
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
##    This function computes the Seasonal Sens's Slope.
##
        na.fail(x)
        if (!is.ts(x)){
            stop("error: input must be ts object")
        }
        fr <- frequency(x)
        if (fr < 2){
            stop("error: minimum of two seasons required.")
        }

        n <- length(x)
        retval <- list(b.sen = NULL,
                   intercept = NULL, nobs = n, method="SeaSLP")
        y <- NULL
        d <- NULL
        varS <- 0
        for (i in 1:fr) {
            y <- x[cycle(x) == i]
            res <- sens.slope(ts(y))
            D.m <- res$D
            d <- c(d, D.m)
            varS <- varS + varS
        }

        b.sen <- median(d, na.rm=TRUE)
        intercepts <- as.vector(x) - b.sen * (1:n)
        retval$intercept <- median(intercepts)
        retval$b.sen <- b.sen
        class(retval) <- "trend.test"
        return(retval)

    }
