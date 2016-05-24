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

########### Country functions ########### 
#' @templateVar name Country
#' @template api
NULL

## GET /admin/countries.json
## Receive a list of all Countries
#' @rdname Country
getCountries <- function(...) {
    .request("countries", ...)$countries
}

## GET /admin/countries/count.json
## Receive a count of all Countries
#' @rdname Country
getCountriesCount <- function(...) {
    .request(.url("countries","count"), ...)$count
}

## GET /admin/countries/#{id}.json
## Receive a single Country
#' @rdname Country
getCountry <- function(countryId, ...) {
    .request(.url("countries",countryId), ...)$country
}

## POST /admin/countries.json
## Create a new Country
#' @rdname Country
createCountry <- function(country, ...) {
    country <- .wrap(country, "country", check=FALSE)
    .request("countries", reqType="POST", data=country, ...)$country
}

## PUT /admin/countries/#{id}.json
## Modify an existing Country
#' @rdname Country
modifyCountry <- function(country, ...) {
    country <- .wrap(country, "country")
    .request(.url("countries",country$country$id), reqType="PUT", data=country, ...)$country
}

## DELETE /admin/countries/#{id}.json
## Remove a Country from the database
#' @rdname Country
deleteCountry <- function(countryId, ...) {
    .request(.url("countries",countryId), reqType="DELETE", ...)
}