#
#   shopifyr: An R Interface to the Shopify API
#
#   Copyright (C) 2014 Charlie Friedemann cfriedem @ gmail.com
#   Shopify API (c) 2006-2014 Shopify Inc.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

########### Refund functions ########### 
#' @param orderId an Order id number
#' @param refundId a Refund id number
#' @templateVar name Refund
#' @templateVar default.params FALSE
#' @templateVar default.return FALSE
#' @return a list corresponding to a Refund
#' @template api
NULL

## GET /admin/orders/#{id}/refunds/#{id}.json
## Receive a single Refund
#' @rdname Refund
getRefund <- function(orderId, refundId, ...) {
    .request(.url("orders",orderId,"refunds",refundId), ...)$refund
}