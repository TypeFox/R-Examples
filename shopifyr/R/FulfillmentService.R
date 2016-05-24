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

########### FulfillmentService functions ########### 
#' @templateVar name FulfillmentService
#' @template api
NULL

## GET /admin/fulfillment_services.json
## Receive a list of all FulfillmentServices
#' @rdname FulfillmentService
getFulfillmentServices <- function(...) {
    .request("fulfillment_services", ...)$fulfillment_services
}

## POST /admin/fulfillment_services.json
## Create a new FulfillmentService
#' @rdname FulfillmentService
createFulfillmentService <- function(fulfillmentService, ...) {
    fulfillmentService <- .wrap(fulfillmentService, "fulfillment_service", check=FALSE)
    .request("fulfillment_services", reqType="POST", data=fulfillmentService, ...)$fulfillment_service
}

## GET /admin/fulfillment_services/#{id}.json
## Receive a single FulfillmentService
#' @rdname FulfillmentService
getFulfillmentService <- function(fulfillmentServiceId, ...) {
    .request(.url("fulfillment_services",fulfillmentServiceId), ...)$fulfillment_service
}

## PUT /admin/fulfillment_services/#{id}.json
## Modify an existing FulfillmentService
#' @rdname FulfillmentService
modifyFulfillmentService <- function(fulfillmentService, ...) {
    fulfillmentService <- .wrap(fulfillmentService, "fulfillment_service")
    .request(.url("fulfillment_services",fulfillmentService$fulfillment_service$id), reqType="PUT", data=fulfillmentService, ...)$fulfillment_service
}

## DELETE /admin/fulfillment_services/#{id}.json
## Remove a FulfillmentService from the database
#' @rdname FulfillmentService
deleteFulfillmentService <- function(fulfillmentServiceId, ...) {
    .request(.url("fulfillment_services",fulfillmentServiceId), reqType="DELETE", ...)
}