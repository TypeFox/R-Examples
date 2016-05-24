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

########### Location functions ########### 
#' @param locationId a Location id number
#' @templateVar name Location
#' @templateVar default.params FALSE
#' @template api
NULL

## GET /admin/locations.json
## Receive a list of all Locations
#' @rdname Location
getLocations <- function(...) {
    .request("locations", ...)$locations
}

## GET /admin/locations/#{id}.json
## Receive a single Location
#' @rdname Location
getLocation <- function(locationId, ...) {
    .request(.url("locations",locationId), ...)$location
}