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

########### Page functions ########### 
#' @templateVar name Page
#' @template api
NULL

## GET /admin/pages.json
## Receive a list of all Pages
#' @rdname Page
getPages <- function(...) {
    .fetchAll("pages", ...)
}

## GET /admin/pages/count.json
## Receive a count of all Pages
#' @rdname Page
getPagesCount <- function(...) {
    .request(.url("pages","count"), ...)$count
}

## GET /admin/pages/#{id}.json
## Receive a single Page
#' @rdname Page
getPage <- function(pageId, ...) {
    .request(.url("pages",pageId), ...)$page
}

## POST /admin/pages.json
## Create a new Page
#' @rdname Page
createPage <- function(page, ...) {
    page <- .wrap(page, "page", check=FALSE)
    .request("pages", reqType="POST", data=page, ...)$page
}

## PUT /admin/pages/#{id}.json
## Modify an existing Page
#' @rdname Page
modifyPage <- function(page, ...) {
    page <- .wrap(page, "page")
    .request(.url("pages",page$page$id), reqType="POST", data=page, ...)$page
}

## DELETE /admin/pages/#{id}.json
## Remove a Page from the database
#' @rdname Page
deletePage <- function(pageId, ...) {
    .request(.url("pages",pageId), reqType="DELETE", ...)
}
