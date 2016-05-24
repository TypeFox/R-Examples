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

########### ProductVariant functions ########### 
#' @param productId a Product id number
#' @templateVar name ProductVariant
#' @templateVar slug variant
#' @templateVar urlSlug product_variant
#' @template api
NULL

## GET /admin/products/#{id}/variants.json
## Receive a list of all Product Variants
#' @rdname ProductVariant
getProductVariants <- function(productId, ...) {
    .fetchAll(.url("products",productId,"variants"), "variants", ...)
}

## GET /admin/products/#{id}/variants/count.json
## Receive a count of all Product Variants
#' @rdname ProductVariant
getProductVariantsCount <- function(productId, ...) {
    .request(.url("products",productId,"variants","count"), ...)$count
}

## GET /admin/variants/#{id}.json
## Receive a single Product Variant
#' @rdname ProductVariant
getProductVariant<- function(variantId, ...) {
    .request(.url("variants",variantId), ...)$variant
}

## POST /admin/products/#{id}/variants.json
## Create a new Product Variant
#' @rdname ProductVariant
createProductVariant<- function(productId, variant, ...) {
    variant <- .wrap(variant, "variant", check=FALSE)
    .request(.url("products",productId,"variants"), reqType="POST", data=variant, ...)$variant
}

## PUT /admin/variants/#{id}.json
## Modify an existing Product Variant
#' @rdname ProductVariant
modifyProductVariant<- function(productId, variant, ...) {
    variant <- .wrap(variant, "variant")
    .request(.url("variants", variant$variant$id), reqType="PUT", data=variant, ...)$variant
}

## DELETE /admin/products/#{id}/variants/#{id}.json
## Remove a Product Variant from the database
#' @rdname ProductVariant
deleteProductVariant<- function(productId, variantId, ...) {
    .request(.url("products",productId,"variants", variantId), reqType="DELETE", ...)
}
