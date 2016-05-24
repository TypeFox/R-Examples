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

########### Event functions ########### 
#' @param eventId a Event id number
#' @templateVar name Event
#' @templateVar default.params FALSE
#' @template api
NULL

## GET /admin/events.json
## Receive a list of all Events
#' @rdname Event
getEvents <- function(...) {
    .fetchAll("events", ...)
}

## GET /admin/events/#{id}.json
## Receive a single Event
#' @rdname Event
getEvent <- function(eventId, ...) {
    .request(.url("events",eventId), ...)$event
}

## GET /admin/events/count.json
## Receive a count of all Events
#' @rdname Event
getEventsCount <- function(...) {
    .request(.url("events","count"), ...)$count
}