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

########### Province functions ########### 
#' @param countryId a Country id number
#' @templateVar name Province
#' @template api
NULL

## GET /admin/countries/#{id}/provinces.json
## Receive a list of all Provinces
#' @rdname Province
getProvinces <- function(countryId, ...) {
    .request(.url("countries",countryId,"provinces"), ...)$provinces
}

## GET /admin/countries/#{id}/provinces/count.json
## Receive a count of all Provinces
#' @rdname Province
getProvincesCount <- function(countryId, ...) {
    .request(.url("countries",countryId,"provinces","count"), ...)$count
}

## GET /admin/countries/#{id}/provinces/#{id}.json
## Receive a single Province
#' @rdname Province
getProvince <- function(countryId, provinceId, ...) {
    .request(.url("countries",countryId,"provinces",provinceId), ...)$province
}

## PUT /admin/countries/#{id}/provinces/#{id}.json
## Modify an existing Province
#' @rdname Province
modifyProvince <- function(countryId, province, ...) {
    province <- .wrap(province, "province")
    .request(.url("countries",countryId,"provinces",province$province$id), reqType="PUT", data=province, ...)$province
}
