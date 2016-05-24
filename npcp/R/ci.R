#################################################################################
##
##   R package npcp by Ivan Kojadinovic Copyright (C) 2014
##
##   This file is part of the R package npcp.
##
##   The R package npcp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package npcp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package npcp. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################


#################################################################################
## Confidence intervals for U-statistics
#################################################################################

## ciU <- function(x, ##method = c("mult", "asym.var"),
##                 statistic = c("variance", "gini", "kendall"),
##                 b = 1, weights = c("parzen", "bartlett"),
##                 N = 1000, init.seq = NULL)
## {
##     ##method <- match.arg(method)
##     statistic <- match.arg(statistic)
##     weights <- match.arg(weights)

##     stopifnot(is.matrix(x))
##     n <- nrow(x)

##     ## kernel
##     h.func <- switch(statistic,
##                      variance = function(x, y) (x - y)^2/2,
##                      gini = function(x, y) abs(x - y),
##                      kendall = function(x, y) prod(x < y) + prod(y < x))

##     #h <- outer(1:n, 1:n, h.func)
##     h <- matrix(0,n,n)
##     for (i in seq_len(n))
##         for (j in seq_len(n))
##             if (i < j)
##             {
##                 h[i,j] <- h.func(x[i,],x[j,])
##                 h[j,i] <- h[i,j]
##             }

##     influ <- colSums(h) / (n-1) ## h1.n without centering term

##     if (is.null(b))
##         b <- bOptU(influ, weights=weights)

##     ## m <- switch(method,
##     ##             "mult" = 1,
##     ##             "asym.var" = 2)

##     ##if (method == "mult")
##     ##{
##         ## initial standard normal sequence for generating dependent multipliers
##     if (is.null(init.seq))
##         init.seq <- rnorm(N * (n + 2 * (b - 1)))
##     else
##         stopifnot(length(init.seq) == N * (n + 2 * (b - 1)))
##     ##}

##     ## test
##     out <- .C("ciU",
##               as.double(h),
##               as.integer(n),
##               as.double(influ),
##               u = double(1),
##               as.integer(N),
##               as.integer(weights == "bartlett"),
##               as.integer(b),
##               u0 = double(N),
##               avar = double(1),
##               avar0 = double(N),
##               as.double(init.seq),
##               PACKAGE = "npcp")

##     list(stat=out$u, stat0=out$u0, b=b, sd=2*sqrt(out$avar), sd0=2*sqrt(out$avar0))
## }
