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

########### OrderRisks functions ########### 
#' @param orderId an Order id number
#' @templateVar name OrderRisks
#' @templateVar slug orderRisk
#' @templateVar urlSlug order_risks
#' @template api
NULL

## POST /admin/orders/#{id}/risks.json
## Create a new Order Risks
#' @rdname OrderRisks
createOrderRisk <- function(orderId, orderRisk, ...) {
    orderRisk <- .wrap(orderRisk, "risk", check=FALSE)
    .request(.url("orders",orderId,"risks"), reqType="POST", data=orderRisk, ...)$risk
}

## GET /admin/orders/#{id}/risks.json
## Receive a list of all Order Risks
#' @rdname OrderRisks
getOrderRisks <- function(orderId, ...) {
    .request(.url("orders",orderId,"risks"), ...)$risks
}

## GET /admin/orders/#{id}/risks/#{id}.json
## Receive a single Order Risk
#' @rdname OrderRisks
getOrderRisk <- function(orderId, orderRiskId, ...) {
    .request(.url("orders",orderId,"risks",orderRiskId), ...)$risk
}

## PUT /admin/orders/#{id}/risks/#{id}.json
## Modify an existing Order Risks
#' @rdname OrderRisks
modifyOrderRisk <- function(orderId, orderRisk, ...) {
    orderRisk <- .wrap(orderRisk, "risk")
    .request(.url("orders",orderId,"risks",orderRisk$risk$id), reqType="PUT", data=orderRisk, ...)$risk
}

## DELETE /admin/orders/#{id}/risks/#{id}.json
## Remove a Order Risks from the database
#' @rdname OrderRisks
deleteOrderRisk <- function(orderId, orderRiskId, ...) {
    .request(.url("orders",orderId,"risks",orderRiskId), reqType="DELETE", ...)
}