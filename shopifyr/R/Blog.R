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

########### Blog functions ########### 
#' @templateVar name Blog
#' @template api
NULL

## GET /admin/blogs.json
## Receive a list of all Blogs
#' @rdname Blog
getBlogs <- function(...) {
    .request("blogs", ...)$blogs
}

## GET /admin/blogs/count.json
## Receive a count of all Blogs
#' @rdname Blog
getBlogsCount <- function(...) {
    .request(.url("blogs","count"), ...)$count
}

## GET /admin/blogs/#{id}.json
## Receive a single Blog
#' @rdname Blog
getBlog <- function(blogId, ...) {
    .request(.url("blogs",blogId), ...)$blog
}

## POST /admin/blogs.json
## Create a new Blog
#' @rdname Blog
createBlog <- function(blog, ...) {
    blog <- .wrap(blog, "blog", check=FALSE)
    .request("blogs", reqType="POST", data=blog, ...)$blog
}

## PUT /admin/blogs/#{id}.json
## Modify an existing Blog
#' @rdname Blog
modifyBlog <- function(blog, ...) {
    blog <- .wrap(blog, "blog")
    .request(.url("blogs",blog$blog$id), reqType="PUT", data=blog, ...)$blog
}

## DELETE /admin/blogs/#{id}.json
## Remove a Blog from the database
#' @rdname Blog
deleteBlog <- function(blogId, ...) {
    .request(.url("blogs",blogId), reqType="DELETE", ...)
}