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

########### Order functions ########### 
#' @templateVar name Order
#' @template api
NULL

## GET /admin/orders.json
## Receive a list of all Orders
#' @rdname Order
getOrders <- function(...) {
    .fetchAll("orders", ...)
}

## GET /admin/orders/#{id}.json
## Receive a single Order
#' @rdname Order
getOrder <- function(orderId, ...) {
    .request(.url("orders",orderId), ...)$order
}

## GET /admin/orders/count.json
## Receive a count of all Orders
#' @rdname Order
getOrdersCount <- function(...) {
    .request(.url("orders","count"), ...)$count
}

## POST /admin/orders/#{id}/close.json
## Close an Order
#' @rdname Order
closeOrder <- function(orderId, ...) {
    .request(.url("orders",orderId,"close"), reqType="POST", data=list(), ...)$order
}

## POST /admin/orders/#{id}/open.json
## Re-open a closed Order
#' @rdname Order
openOrder <- function(orderId, ...) {
    .request(.url("orders",orderId,"open"), reqType="POST", data=list(), ...)$order
}

## POST /admin/orders/#{id}/cancel.json
## Cancel an Order
#' @rdname Order
cancelOrder <- function(orderId, ...) {
    .request(.url("orders",orderId,"cancel"), reqType="POST", data=list(), ...)$order
}

## POST /admin/orders.json
## Create a new Order
#' @rdname Order
createOrder <- function(order, ...) {
    order <- .wrap(order, "order", check=FALSE)
    .request("orders", reqType="POST", data=order, ...)$order
}

## PUT /admin/orders/#{id}.json
## Modify an existing Order
#' @rdname Order
modifyOrder <- function(order, ...) {
    order <- .wrap(order, "order")
    .request(.url("orders",order$order$id), reqType="PUT", data=order, ...)$order
}

## DELETE /admin/orders/#{id}.json
## Remove a Order from the database
#' @rdname Order
deleteOrder <- function(orderId, ...) {
    .request(.url("orders",orderId), reqType="DELETE", ...)
}