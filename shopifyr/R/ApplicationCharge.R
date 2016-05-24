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

########### ApplicationCharge functions ########### 
#' @templateVar slug charge
#' @templateVar name ApplicationCharge
#' @template api
NULL

## POST /admin/application_charges.json
## Create a new one-time application charge
#' @rdname ApplicationCharge
createApplicationCharge <- function(charge, ...) {
    charge <- .wrap(charge, "application_charge", check=FALSE)
    .request("application_charges", reqType="POST", data=charge, ...)$application_charge
}

## GET /admin/application_charges/#{id}.json
## Receive a single ApplicationCharge
#' @rdname ApplicationCharge
getApplicationCharge <- function(chargeId, ...) {
    .request(.url("application_charges",chargeId), ...)$application_charge
}

## GET /admin/application_charges.json
## Retrieve all one-time application charges
#' @rdname ApplicationCharge
getApplicationCharges <- function(...) {
    .request("application_charges", ...)$application_charges
}

## POST /admin/application_charges/#{id}/activate.json
## Activate a one-time application charge
#' @rdname ApplicationCharge
activateApplicationCharge <- function(charge, ...) {
    charge <- .wrap(charge, "application_charge")
    .request(.url("application_charges",charge$application_charge$id), reqType="POST", data=charge, ...)$application_charge
}