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

### This demo shows how the HAC package can be used for estimationg NACs

require(HAC)
require(copula)

## Build an 'nacopula' object (nested Archimedean copula (NAC))
theta <- 2:5
copG <- onacopulaL("Gumbel", list(theta[1], NULL, list(list(theta[2], c(2,1)),
                                                       list(theta[4], c(5,6)),
                                                       list(theta[3], c(4,3)))))
## Sample from copG
set.seed(271)
U <- rnacopula(1000, copula=copG)

## fitCopula(copG, U) does not provide fitting capabilities for HACs/NACs yet
## but we can convert copG to a 'hac' object
hacG <- nacopula2hac(copG)
plot(hacG) # plot method

## Parameters can either be estimated based on a fixed structure...
colnames(U) <- paste(1:ncol(U))
hac.fixed <- estimate.copula(U, hac=hacG)
## ... or the structure can be estimated as well:
hac.flex <- estimate.copula(U, type=hacG$type)

## Show the estimates
plot(hac.fixed)
plot(hac.flex)

## Last but not least, the estimation results can be re-converted
## into 'nacopula'-objects again
cop.fixed <- hac2nacopula(hac.fixed)
cop.flex <- hac2nacopula(hac.flex)
