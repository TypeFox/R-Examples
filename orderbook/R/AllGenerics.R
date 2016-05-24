################################################################################
##
##
## orderbook.R: Generic functions for the limitob class
##
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
## along with limitob.  If not, see <http://www.gnu.org/licenses/>.

################################################################################

.onLoad <- function(lib, pkg) require(methods)

setClass("orderbook", representation(current.ob   = "data.frame",
                                     current.time = "numeric",
                                     feed         = "character",
                                     feed.index   = "numeric",
                                     ob.data      = "hash",
                                     trade.data   = "vector",
                                     trader       = "logical"
                                     ))

if(!isGeneric("spread"))
    setGeneric("spread", function(object, ...) standardGeneric("spread"))

if(!isGeneric("display"))
    setGeneric("display", function(object, n = 5, short = TRUE, ...)
               standardGeneric("display"))

if(!isGeneric("ask.price.levels"))
    setGeneric("ask.price.levels", function(object, ...)
               standardGeneric("ask.price.levels"))

if(!isGeneric("bid.price.levels"))
    setGeneric("bid.price.levels", function(object, ...)
               standardGeneric("bid.price.levels"))

if(!isGeneric("total.price.levels"))
    setGeneric("total.price.levels", function(object, ...)
               standardGeneric("total.price.levels"))

if(!isGeneric("bid.orders"))
    setGeneric("bid.orders", function(object, ...)
               standardGeneric("bid.orders"))

if(!isGeneric("ask.orders"))
    setGeneric("ask.orders", function(object, ...)
               standardGeneric("ask.orders"))

if(!isGeneric("total.orders"))
    setGeneric("total.orders", function(object, ...)
               standardGeneric("total.orders"))

if(!isGeneric("mid.point"))
    setGeneric("mid.point", function(object, ...)
               standardGeneric("mid.point"))

if(!isGeneric("inside.market"))
    setGeneric("inside.market", function(object, invis = FALSE, ...)
               standardGeneric("inside.market"))

if(!isGeneric("best.bid"))
    setGeneric("best.bid", function(object, ...)
			   standardGeneric("best.bid"))

if(!isGeneric("best.ask"))
    setGeneric("best.ask", function(object, ...)
	           standardGeneric("best.ask"))

if(!isGeneric("add.order"))
    setGeneric("add.order", function(object, price, size, type,
               time=NULL, id = NULL, status = FALSE)
               standardGeneric("add.order"))

if(!isGeneric("remove.order"))
    setGeneric("remove.order", function(object, id)
               standardGeneric("remove.order"))

if(!isGeneric("replace.order"))
    setGeneric("replace.order", function(object, id, size)
               standardGeneric("replace.order"))

if (!isGeneric("market.order"))
    setGeneric("market.order", function(object, size, type)
               standardGeneric("market.order"))

if(!isGeneric("reset"))
    setGeneric("reset", function(object)
               standardGeneric("reset"))

if(!isGeneric("read.orders"))
    setGeneric("read.orders", function(object, n = 1000)
               standardGeneric("read.orders"))

if(!isGeneric("read.time"))
    setGeneric("read.time", function(object, n)
               standardGeneric("read.time"))

if(!isGeneric("next.trade"))
    setGeneric("next.trade", function(object)
               standardGeneric("next.trade"))

if(!isGeneric("previous.trade"))
    setGeneric("previous.trade", function(object)
               standardGeneric("previous.trade"))

if(!isGeneric("midpoint.return"))
    setGeneric("midpoint.return", function(object, tradenum, time)
               standardGeneric("midpoint.return"))

if(!isGeneric("view.trade"))
    setGeneric("view.trade", function(object, tradenum)
               standardGeneric("view.trade"))
