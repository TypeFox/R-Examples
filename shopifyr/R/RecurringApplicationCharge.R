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

########### RecurringApplicationCharge functions ########### 
#' @templateVar name RecurringApplicationCharge
#' @templateVar slug charge
#' @template api
NULL

## POST /admin/recurring_application_charges.json
## Create a recurring application charge
#' @rdname RecurringApplicationCharge
createRecurringCharge <- function(charge, ...) {
    charge <- .wrap(charge, "recurring_application_charge", check=c("price","name"))
    .request("recurring_application_charges", reqType="POST", data=charge, ...)$recurring_application_charge
}

## GET /admin/recurring_application_charges/#{id}.json
## Receive a single RecurringApplicationCharge
#' @rdname RecurringApplicationCharge
getRecurringCharge <- function(chargeId, ...) {
    .request(.url("recurring_application_charges",chargeId), ...)$recurring_application_charge
}

## GET /admin/recurring_application_charges.json
## Retrieve all recurring application charges
#' @rdname RecurringApplicationCharge
getRecurringCharges <- function(...) {
    .request("recurring_application_charges", ...)$recurring_application_charges
}

## POST /admin/recurring_application_charges/#{id}/activate.json
## Activate a recurring application charge
#' @rdname RecurringApplicationCharge
activateRecurringCharge <- function(charge, ...) {
    charge <- .wrap(charge, "recurring_application_charge")
    .request(.url("recurring_application_charges",charge$charge$id,"activate"), reqType="POST", data=charge, ...)
}

## DELETE /admin/recurring_application_charges/#{id}.json
## Cancel a recurring application charge
#' @rdname RecurringApplicationCharge
cancelRecurringCharge <- function(chargeId, ...) {
    .request(.url("recurring_application_charges",chargeId), reqType="DELETE", ...)
}