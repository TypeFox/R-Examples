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

########### Webhook functions ########### 
#' @templateVar name Webhook
#' @template api
NULL

## GET /admin/webhooks.json
## Receive a list of all Webhooks
#' @rdname Webhook
getWebhooks <- function(...) {
    .fetchAll("webhooks", ...)
}

## GET /admin/webhooks/count.json
## Receive a count of all Webhooks
#' @rdname Webhook
getWebhooksCount <- function(...) {
    .request(.url("webhooks","count"), ...)$count
}

## GET /admin/webhooks/#{id}.json
## Receive a single Webhook
#' @rdname Webhook
getWebhook <- function(webhookId, ...) {
    .request(.url("webhooks",webhookId), ...)$webhook
}

## POST /admin/webhooks.json
## Create a new Webhook
#' @rdname Webhook
createWebhook <- function(webhook, ...) {
    webhook <- .wrap(webhook, "webhook", check=c("address","topic"))
    .request("webhooks", reqType="POST", data=webhook, ...)$webhook
}

## PUT /admin/webhooks/#{id}.json
## Modify an existing Webhook
#' @rdname Webhook
modifyWebhook <- function(webhook, ...) {
    webhook <- .wrap(webhook, "webhook")
    .request(.url("webhooks",webhook$webhook$id), reqType="PUT", data=webhook, ...)$webhook
}

## DELETE /admin/webhooks/#{id}.json
## Remove a Webhook from the database
#' @rdname Webhook
deleteWebhook <- function(webhookId, ...) {
    .request(.url("webhooks",webhookId), reqType="DELETE", ...)
}