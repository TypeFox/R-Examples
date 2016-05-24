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

########### CustomCollection functions ########### 
#' @templateVar name CustomCollection
#' @templateVar urlSlug custom_collection
#' @template api
NULL

## GET /admin/custom_collections.json
## Receive a list of all CustomCollections
#' @rdname CustomCollection
getCustomCollections <- function(...) {
    .fetchAll("custom_collections", ...)
}

## GET /admin/custom_collections/count.json
## Receive a count of all CustomCollections
#' @rdname CustomCollection
getCustomCollectionsCount <- function(...) {
    .request(.url("custom_collections","count"), ...)$count
}

## GET /admin/custom_collections/#{id}.json
## Receive a single CustomCollection
#' @rdname CustomCollection
getCustomCollection <- function(customCollectionId, ...) {
    .request(.url("custom_collections",customCollectionId), ...)$custom_collection
}

## POST /admin/custom_collections.json
## Create a new CustomCollection
#' @rdname CustomCollection
createCustomCollection <- function(customCollection, ...) {
    customCollection <- .wrap(customCollection, "custom_collection", check=FALSE)
    .request("custom_collections", reqType="POST", data=customCollection, ...)$custom_collection
}

## PUT /admin/custom_collections/#{id}.json
## Modify an existing CustomCollection
#' @rdname CustomCollection
modifyCustomCollection <- function(customCollection, ...) {
    customCollection <- .wrap(customCollection, "custom_collection")
    .request(.url("custom_collections",customCollection$custom_collection$id), reqType="PUT", data=customCollection, ...)$custom_collection
}

## DELETE /admin/custom_collections/#{id}.json
## Remove a CustomCollection from the database
#' @rdname CustomCollection
deleteCustomCollection <- function(customCollectionId, ...) {
    .request(.url("custom_collections",customCollectionId), reqType="DELETE", ...)
}
