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

########### Asset functions ########### 
#' @param themeId a Theme id number
#' @param assetKey an Asset key e.g. \code{"templates/index.liquid"}
#' @param asset a  list containing Asset fields
#' @templateVar name Asset
#' @templateVar article an
#' @templateVar default.params FALSE
#' @template api
NULL

## GET /admin/themes/#{id}/assets.json
## Receive a list of all Assets
#' @rdname Asset
getAssets <- function(themeId, ...) {
    .request(.url("themes",themeId,"assets"), ...)$assets
}

## GET /admin/themes/#{id}/assets.json?asset[key]=templates/index.liquid&theme_id=828155753
## Receive a single Asset
#' @rdname Asset
getAsset <- function(themeId, assetKey, ...) {
    .request(.url("themes",themeId,"assets"), `asset[key]`=assetKey, theme_id=themeId, ...)$asset
}

## PUT /admin/themes/#{id}/assets.json
## Creating or Modifying an Asset
#' @rdname Asset
#' @aliases modifyAsset
createAsset <- modifyAsset <- function(themeId, asset, ...) {
    asset <- .wrap(asset, "asset", check="key")
    .request(.url("themes",themeId,"assets"), reqType="PUT", data=asset, ...)$asset
}

## DELETE /admin/themes/#{id}/assets.json?asset[key]=assets/bg-body.gif
## Remove a Asset from the database
#' @rdname Asset
deleteAsset <- function(themeId, assetKey, ...) {
    .request(.url("themes",themeId,"assets"), `asset[key]`=assetKey, theme_id=themeId, reqType="DELETE", ...)
}