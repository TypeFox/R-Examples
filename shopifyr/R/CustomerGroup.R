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

########### CustomerGroup functions ########### 
#' @templateVar name CustomerGroup
#' @template api
#' @aliases CustomerSavedSearch
NULL

## GET /admin/customer_saved_searches.json
## Receive a list of all CustomerGroups
#' @rdname CustomerGroup
getCustomerGroups <- function(...) {
    .fetchAll("customer_saved_searches", ...)
}

## GET /admin/customer_saved_searches/count.json
## Receive a count of all CustomerGroups
#' @rdname CustomerGroup
getCustomerGroupsCount <- function(...) {
    .request(.url("customer_saved_searches","count"), ...)$count
}

## GET /admin/customer_saved_searches/#{id}.json
## Receive a single CustomerGroup
#' @rdname CustomerGroup
getCustomerGroup <- function(customerGroupId, ...) {
    .request(.url("customer_saved_searches",customerGroupId), ...)$customer_saved_search
}

## GET /admin/customer_saved_searches/#{id}/customers.json
## Receive all Customers belonging to a CustomerGroup
#' @rdname CustomerGroup
getCustomerGroupCustomers <- function(customerGroupId, ...) {
    .fetchAll(.url("customer_saved_searches",customerGroupId,"customers"), "customers", ...)
}

## POST /admin/customer_saved_searches.json
## Create a new CustomerGroup
#' @rdname CustomerGroup
createCustomerGroup <- function(customerGroup, ...) {
    customerGroup <- .wrap(customerGroup, "customer_saved_search", check=FALSE)
    .request("customer_saved_searches", reqType="POST", data=customerGroup, ...)$customer_saved_search
}

## PUT /admin/customer_saved_searches/#{id}.json
## Modify an existing CustomerGroup
#' @rdname CustomerGroup
modifyCustomerGroup <- function(customerGroup, ...) {
    customerGroup <- .wrap(customerGroup, "customer_saved_search")
    .request(.url("customer_saved_searches",customerGroup$customer_saved_searches$id), reqType="PUT", data=customerGroup, ...)$customer_saved_search
}

## DELETE /admin/customer_saved_searches/#{id}.json
## Remove a CustomerGroup from the database
#' @rdname CustomerGroup
deleteCustomerGroup <- function(customerGroupId, ...) {
    .request(.url("customer_saved_search",customerGroupId), reqType="DELETE", ...)
}