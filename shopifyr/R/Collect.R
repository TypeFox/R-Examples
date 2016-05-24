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

########### Collect functions ########### 
#' @templateVar name Collect
#' @template api
NULL

## POST /admin/collects.json
## Create a new Collect
#' @rdname Collect
createCollect <- function(collect, ...) {
    collect <- .wrap(collect, "collect", check=c("product_id","collection_id"))
    .request("collects", reqType="POST", data=collect, ...)$collect
}

## DELETE /admin/collects/#{id}.json
## Remove a Collect from the database
#' @rdname Collect
deleteCollect <- function(collectId, ...) {
    .request(.url("collects",collectId), reqType="DELETE", ...)
}

## GET /admin/collects.json
## Receive a list of all Collects
#' @rdname Collect
getCollects <- function(...) {
    .fetchAll("collects", ...)
}

## GET /admin/collects/count.json
## Receive a count of all Collects
#' @rdname Collect
getCollectsCount <- function(...) {
    .request(.url("collects","count"), ...)$count
}

## GET /admin/collects/#{id}.json
## Receive a single Collect
#' @rdname Collect
getCollect <- function(collectId, ...) {
    .request(.url("collects",collectId), ...)$collect
}