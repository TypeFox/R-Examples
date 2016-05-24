## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

### Kemp 1981, Applied Statistics 30(3), pp349--253 for RNG
### The cf is Log[1 - a exp(i t)] / Log[1 - a]
## dlogseries <- function(x, alpha, log = FALSE) {
##   val <-  - alpha^x / x / log(1 - alpha)
##   if (log) log(val) else val
## }

## plogseries <- function(q, alpha, lower.tail = TRUE, log.p = FALSE) {
##   if (!lower.tail) q <- 1 - q
##   val <- 1 - alpha^(1 + q) * hyperg2F1(1 + q, 1, 2 + q, alpha) / (1 + q)
##   if (log.p) log(val) else val
## }

##' only used in  rfrankCopula() currently:
rlogseries <- function(n, alpha) {
  val <- integer(n)
  alpha <- rep(alpha, len = n)
  .C(rlogseries_R, as.integer(n), as.double(alpha), val = as.integer(val))$val
}
