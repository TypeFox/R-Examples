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

########### Redirect functions ########### 
#' @templateVar name Redirect
#' @template api
NULL

## GET /admin/redirects.json
## Receive a list of all Redirects
#' @rdname Redirect
getRedirects <- function(...) {
    .fetchAll("redirects", ...)
}

## GET /admin/redirects/count.json
## Receive a count of all Redirects
#' @rdname Redirect
getRedirectsCount <- function(...) {
    .request(.url("redirects","count"), ...)$count
}

## GET /admin/redirects/#{id}.json
## Receive a single Redirect
#' @rdname Redirect
getRedirect <- function(redirectId, ...) {
    .request(.url("redirects",redirectId), ...)$redirect
}

## POST /admin/redirects.json
## Create a new Redirect
#' @rdname Redirect
createRedirect <- function(redirect, ...) {
    redirect <- .wrap(redirect, "redirect", check=c("path","target"))
    .request("redirects", reqType="POST", data=redirect, ...)$redirect
}

## PUT /admin/redirects/#{id}.json
## Modify an existing Redirect
#' @rdname Redirect
modifyRedirect <- function(redirect, ...) {
    redirect <- .wrap(redirect, "redirect")
    .request(.url("redirects",redirect$redirect$id), reqType="PUT", data=redirect, ...)$redirect
}

## DELETE /admin/redirects/#{id}.json
## Remove a Redirect from the database
#' @rdname Redirect
deleteRedirect <- function(redirectId, ...) {
    .request(.url("redirects",redirectId), reqType="DELETE", ...)
}