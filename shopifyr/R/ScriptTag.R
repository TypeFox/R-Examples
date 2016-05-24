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

########### ScriptTag functions ########### 
#' @templateVar name ScriptTag
#' @template api
NULL

## GET /admin/script_tags.json
## Receive a list of all ScriptTags
#' @rdname ScriptTag
getScriptTags <- function(...) {
    .fetchAll("script_tags", ...)
}

## GET /admin/script_tags/count.json
## Receive a count of all ScriptTags
#' @rdname ScriptTag
getScriptTagsCount <- function(...) {
    .request(.url("script_tags","count"), ...)$count
}

## GET /admin/script_tags/#{id}.json
## Receive a single ScriptTag
#' @rdname ScriptTag
getScriptTag <- function(scriptTagId, ...) {
    .request(.url("script_tags",scriptTagId), ...)$script_tag
}

## POST /admin/script_tags.json
## Create a new ScriptTag
#' @rdname ScriptTag
createScriptTag <- function(scriptTag, ...) {
    scriptTag <- .wrap(scriptTag, "script_tag", check=c("src","event"))
    .request("script_tags", reqType="POST", data=scriptTag, ...)$script_tag
}

## PUT /admin/script_tags/#{id}.json
## Modify an existing ScriptTag
#' @rdname ScriptTag
modifyScriptTag <- function(scriptTag, ...) {
    scriptTag <- .wrap(scriptTag, "script_tag")
    .request(.url("script_tags",scriptTag$script_tag$id), reqType="PUT", data=scriptTag, ...)$script_tag
}

## DELETE /admin/script_tags/#{id}.json
## Remove a ScriptTag from the database
#' @rdname ScriptTag
deleteScriptTag <- function(scriptTagId, ...) {
    .request(.url("script_tags",scriptTagId), reqType="DELETE", ...)
}