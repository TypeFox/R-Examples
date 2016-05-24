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

########### Fulfillment functions ########### 
#' @param orderId an Order id number
#' @templateVar name Fulfillment
#' @template api
NULL

## GET /admin/orders/#{id}/fulfillments.json
## Receive a list of all Fulfillments
#' @rdname Fulfillment
getFulfillments <- function(orderId, ...) {
    .fetchAll(.url("orders",orderId,"fulfillments"), "fulfillments", ...)
}

## GET /admin/orders/#{id}/fulfillments/count.json
## Receive a count of all Fulfillments
#' @rdname Fulfillment
getFulfillmentsCount <- function(orderId, ...) {
    .request(.url("orders",orderId,"fulfillments","count"), ...)$count
}

## GET /admin/orders/#{id}/fulfillments/#{id}.json
## Receive a single Fulfillment
#' @rdname Fulfillment
getFulfillment <- function(orderId, fulfillmentId, ...) {
    .request(.url("orders",orderId,"fulfillments",fulfillmentId), ...)$fulfillment
}

## POST /admin/orders/#{id}/fulfillments.json
## Create a new Fulfillment
#' @rdname Fulfillment
createFulfillment <- function(orderId, fulfillment, ...) {
    fulfillment <- .wrap(fulfillment, "fulfillment", check=FALSE)
    .request(.url("orders",orderId,"fulfillments"), reqType="POST", data=fulfillment, ...)$fulfillment
}

## PUT /admin/orders/#{id}/fulfillments/#{id}.json
## Modify an existing Fulfillment
#' @rdname Fulfillment
modifyFulfillment <- function(orderId, fulfillment, ...) {
    fulfillment <- .wrap(fulfillment, "fulfillment")
    .request(.url("orders",orderId,"fulfillments",fulfillment$fulfillment$id), reqType="PUT", data=fulfillment, ...)$fulfillment
}

## POST /admin/orders/#{id}/fulfillments/#{id}/complete.json
## Complete a pending fulfillment
#' @rdname Fulfillment
completeFulfillment <- function(orderId, fulfillmentId, ...) {
    .request(.url("orders",orderId,"fulfillments",fulfillmentId,"complete"), reqType="POST", data=list(), ...)$fulfillment
}

## POST /admin/orders/#{id}/fulfillments/#{id}/cancel.json
## Cancel a pending fulfillment
#' @rdname Fulfillment
cancelFulfillment <- function(orderId, fulfillmentId, ...) {
    .request(.url("orders",orderId,"fulfillments",fulfillmentId,"cancel"), reqType="POST", data=list(), ...)$fulfillment
}