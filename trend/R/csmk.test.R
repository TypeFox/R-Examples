csmk.test <-
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
##    This function computes the Correlated Seasonal Mann Kendall test.
##    It calls multivar.MK.test()
##
    if (!is.ts(x))
        stop("error: input must be ts object")
    fr <- frequency(x)
    if(fr < 2)
        stop("error: minimum of two seasons required.")
    ts.tsp <- tsp(x)
    n <- length(abs(ts.tsp[1]):abs(ts.tsp[2]))
    y <- matrix(data=NA, nrow=n, ncol=fr)
    for (i in 1:fr) {
        y[,i] <- x[cycle(x) == i]
    }
    res <- multivar.MK.test(y, method="CSMK")
    class(res) <- "trend.test"
    res
}

