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

########### Metafield functions ########### 
#' @param resourceName the name of a resource e.g. \code{"products", "smart_collections", "product_image"}
#' @param resourceId the id number of the resource, if applicable (for example Shop has no id)
#' @templateVar name Metafield
#' @template api
NULL

## GET /admin/metafields.json
## Get metafields that belong to a store
## GET /admin/metafields.json?metafield[owner_id]=850703190&metafield[owner_resource]=product_image
## Get metafields that belong to a product image
#' @rdname Metafield
getMetafields <- function(resourceName, resourceId = NULL, ...) {
    if (resourceName == "shop") resourceName <- NULL
    #.fetchAll(.url(resourceName,resourceId,"metafields"), "metafields", ...)$metafield # doesn't work for all resources
    .fetchAll("metafields", `metafield[owner_resource]`=resourceName, `metafield[owner_id]`=resourceId, ...)
}

## GET /admin/metafields/count.json
## Get a count of metafields that belong to a store
## GET /admin/products/#{id}/metafields/count.json
## Get a count of metafields that belong to a product
#' @rdname Metafield
getMetafieldsCount <- function(resourceName, resourceId = NULL, ...) {
    #.request(.url(resourceName,resourceId,"metafields","count"), ...) # doesn't work for all resources
    .request(.url("metafields","count"), `metafield[owner_resource]`=resourceName, `metafield[owner_id]`=resourceId, ...)$count
}

## GET /admin/metafields/#{id}.json
## Get a single store metafield by its ID
## GET /admin/products/#{id}/metafields/#{id}.json
## Get a single product metafield by its ID
#' @rdname Metafield
getMetafield <- function(metafieldId, ...) {
    .request(.url("metafields",metafieldId), ...)$metafield
}

## POST /admin/metafields.json
## Create a new metafield for a store
## POST /admin/products/#{id}/metafields.json
## Create a new metafield for a product
#' @rdname Metafield
createMetafield <- function(resourceName, resourceId = NULL, metafield, ...) {
    metafield <- .wrap(metafield, "metafield", check="key")
    if (resourceName == "shop") resourceName <- NULL
    .request(.url(resourceName,resourceId,"metafields"), reqType="POST", data=metafield, ...)$metafield
}

## PUT /admin/metafields/#{id}.json
## Update a store metafield
## PUT /admin/products/#{id}/metafields/#{id}.json
## Update a product metafield
#' @rdname Metafield
modifyMetafield <- function(metafield, ...) {
    metafield <- .wrap(metafield, "metafield")
    .request(.url("metafields",metafield$metafield$id), reqType="PUT", data=metafield, ...)$metafield
}

## DELETE /admin/metafields/#{id}.json
## Delete a store metafield
## DELETE /admin/products/#{id}/metafields/#{id}.json
## Delete a product metafield
#' @rdname Metafield
deleteMetafield <- function(metafieldId, ...) {
    .request(.url("metafields",metafieldId), reqType="DELETE", ...)
}