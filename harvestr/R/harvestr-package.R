{###############################################################################
# harvestr-packages.R
# This file is part of the R package harvestr.
# 
# Copyright 2012 Andrew Redd
# Date: 6/2/2012
# 
# DESCRIPTION
# ===========
# package documentation
# 
# LICENSE
# ========
# harvestr is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# dostats. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################


#' \code{harvestr} package
#' @name harvestr
#' @aliases package-harvestr harvestr
#' @docType package
#' @aliases harvestr package-harvestr
#' @title A Simple Reproducible Parallel Simulation Framework
#' @author Andrew Redd <amredd_at_gmail.com>
#' 
#' The harvestr package is a framework for parallel reproducible simulations.
#' 
#' The functions to know about are:
#' \enumerate{
#'   \item \code{\link{gather}} - which gathers parallel seeds.
#'   \item \code{\link{farm}}  - which uses the saved seeds from gather to replicate an expression,
#'         once for each seed.
#'   \item \code{\link{harvest}} -  which uses objects from farm, that have saved seed attributes, 
#'         to continue evaluation from where farm finished.
#'   \item \code{\link{reap}} - is used by harvest for a single item
#'   \item \code{\link{plant}} - is used to set seeds for a list of predefined objects so that harvest
#'         can be used on it.
#'   \item \code{\link{sprout}} - Generate independent sub-streams.
#'   \item \code{\link{graft}} - Replicate and object in independent substreams
#'                               of random numbers.
#' }
#' 
#' @section Caching:
#' The functions in \code{harvestr} can cache results for faster and
#' interuptible simulations.  This option defaults to \code{FALSE} but can be 
#' chosen by specifying the \code{cache} parameter in any of the functions
#' that produce results.
#' 
#' The caching is performed by saving a RData file in a specified caching
#' directory.  The default directory is named "harvestr-cache" and resides 
#' under the \link[base:getwd]{working directory}.  This can be specified by setting
#' the \code{harvestr.cache.dir} \code{\link{option}}.  Files in this directory 
#' use file names derived from hashes of the expression to evaluate.  Do not
#' modify the file names.
#' 
#' @section Options:
#' The following options control behavior and default values for harvestr.
#' \enumerate{
#'   \item \code{harvestr.use.cache=FALSE} - Should results be cached for fault 
#'         tollerance and accelerated reproducability?
#'   \item \code{harvestr.cache.dir="harvestr-cache"} - The directory to use for
#'          storing cached results.
#'   \item \code{harvestr.time=FALSE} - Should results be timed?
#'   \item \code{harvestr.use.try=TRUE} - Should the vectorized calls use try
#'         to increase fault tollerance?
#'   \item \code{harvestr.try.silent=FALSE} - Should try be run silently?
#'   \item \code{harvestr.try.summary=TRUE} - Print a result if errors were found?
#'   \item \code{harvestr.parallel=FALSE} - Run results in parallel?
#' }
NULL