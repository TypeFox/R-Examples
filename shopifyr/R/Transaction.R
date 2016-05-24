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

########### Transaction functions ########### 
#' @param orderId an Order id number
#' @templateVar name Transaction
#' @template api
NULL

## GET /admin/orders/#{id}/transactions.json
## Receive a list of all Transactions
#' @rdname Transaction
getTransactions <- function(orderId, ...) {
    .request(.url("orders",orderId,"transactions"), ...)$transactions
}

## GET /admin/orders/#{id}/transactions/count.json
## Receive a count of all Transactions
#' @rdname Transaction
getTransactionsCount <- function(orderId, ...) {
    .request(.url("orders",orderId,"transactions","count"), ...)$count
}

## GET /admin/orders/#{id}/transactions/#{id}.json
## Receive a single Transaction
#' @rdname Transaction
getTransaction <- function(orderId, transactionId, ...) {
    .request(.url("orders",orderId,"transactions",transactionId), ...)$transaction
}

## POST /admin/orders/#{id}/transactions.json
## Create a new Transaction
#' @rdname Transaction
createTransaction <- function(orderId, transaction, ...) {
    transaction <- .wrap(transaction, "transaction", "kind")
    .request(.url("orders",orderId,"transactions"), reqType="POST", data=transaction, ...)$transaction
}