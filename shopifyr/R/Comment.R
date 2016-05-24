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

########### Comment functions ########### 
#' @param blogId a Blog id number (leave blank to fetch all comments for the shop)
#' @param articleId an Article id number (leave blank to fetch all comments for this blog)
#' @templateVar name Comment
#' @template api
NULL

## GET /admin/comments.json?article_id=134645308&blog_id=241253187
## Receive a list of all Comments
#' @rdname Comment
getComments <- function(blogId, articleId, ...) {
    if (missing(blogId)) blogId <- NULL
    if (missing(articleId)) articleId <- NULL
    
    .fetchAll("comments", blog_id=blogId, article_id=articleId, ...)
}

## GET /admin/comments/count.json?article_id=134645308&blog_id=241253187
## Receive a count of all Comments
#' @rdname Comment
getCommentsCount <- function(blogId, articleId, ...) {
    if (missing(blogId)) blogId <- NULL
    if (missing(articleId)) articleId <- NULL
    
    .request(.url("comments","count"), blog_id=blogId, article_id=articleId, ...)$count
}
## GET /admin/comments/#{id}.json
## Receive a single Comment
#' @rdname Comment
getComment <- function(commentId, ...) {
    .request(.url("comments",commentId), ...)$comment
}

## POST /admin/comments.json
## Create a new Comment
#' @rdname Comment
createComment <- function(comment, ...) {
    comment <- .wrap(comment, "comment", check=FALSE)
    .request("comments", reqType="POST", data=comment, ...)$comment
}

## PUT /admin/comments/#{id}.json
## Modify an existing Comment
#' @rdname Comment
modifyComment <- function(comment, ...) {
    comment <- .wrap(comment, "comment")
    .request(.url("comments",comment$comment$id), reqType="PUT", data=comment, ...)$comment
}

## POST /admin/comments/#{id}/spam.json
## Mark a Comment as spam
#' @rdname Comment
markCommentAsSpam <- function(commentId, ...) {
    .request(.url("comments",commentId,"spam"), reqType="POST", data=list())$comment
}

## POST /admin/comments/#{id}/not_spam.json
## Mark a Comment as not spam
#' @rdname Comment
markCommentAsNotSpam <- function(commentId, ...) {
    .request(.url("comments",commentId,"not_spam"), reqType="POST", data=list())$comment
}

## POST /admin/comments/#{id}/approve.json
## Approve a Comment
#' @rdname Comment
approveComment <- function(commentId, ...) {
    .request(.url("comments",commentId,"approve"), reqType="POST", data=list())$comment
}

## POST /admin/comments/#{id}/remove.json
## Remove a Comment
#' @rdname Comment
removeComment <- function(commentId, ...) {
    .request(.url("comments",commentId,"remove"), reqType="POST", data=list())$comment
}

## POST /admin/comments/#{id}/restore.json
## Restore a Comment
#' @rdname Comment
restoreComment <- function(commentId, ...) {
    .request(.url("comments",commentId,"restore"), reqType="POST", data=list())$comment
}