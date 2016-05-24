## read.orders.time.test.R: Tests for read order functions
##
## limitob is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## limitob is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.


library(orderbook)
load("read.orders.time.test.RData")

filename <- system.file("extdata", "sample.txt",
                        package = "orderbook")

ob <- orderbook(file = filename)

# read 5000 orders
ob <- read.orders(ob, 5000)

# 64bit machines will have the data frame in a different order than
# 32bit machines, so a direct call to identical won't work. This
# doesn't matter for any of the functions, so is not a problem.

test <- test[order(test$id),]
rownames(test) <- NULL

temp <- ob@current.ob[order(ob@current.ob$id),]
rownames(temp) <- NULL

stopifnot(isTRUE(identical(test, temp)))


# go back 5 orders and then forward 5 orders
ob <- read.orders(ob, -5)
ob <- read.orders(ob, 5)

test <- test[order(test$id),]
rownames(test) <- NULL

temp <- ob@current.ob[order(ob@current.ob$id),]
rownames(temp) <- NULL

stopifnot(isTRUE(identical(test, temp)))

## Read time, should take us to 4981 orders

ob <- read.time(ob, "9:32:11")
ob <- read.orders(ob, 18)

test <- test[order(test$id),]
rownames(test) <- NULL

temp <- ob@current.ob[order(ob@current.ob$id),]
rownames(temp) <- NULL

stopifnot(isTRUE(identical(test, temp)))

