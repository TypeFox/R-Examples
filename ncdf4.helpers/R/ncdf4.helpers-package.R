#' ncdf4.helpers: helper functions for NetCDF files.
#' 
#' This package provides a number of helper functions for NetCDF files opened using the \code{ncdf4} package.
#' 
#' Dealing with NetCDF format data is unnecessarily difficult. The \code{ncdf4} package does a good job of making many lower-level operations easier. The \code{ncdf4.helpers} package aims to build higher-level functions upon the foundation of \code{ncdf4}.
#' 
#' One concept central to much of the package is the idea of indexing, and dealing with data, by axis rather than by indices or by specific dimension names. The axes used are:
#' \itemize{
#' \item X (the horizontal axis)
#' \item Y (the vertical axis)
#' \item Z (the pressure / depth axis)
#' \item S (the reduced spatial grid axis)
#' \item T (the time axis)
#' }
#' 
#' Indexing by axis avoids the pitfalls of using data in forms other than (X, Y, Z, T), such as (T, X, Y). Avoiding using dimension names directly avoids problems when using projected data. 
#'
#' The functions in the package can be broken down into the following categories:
#' \enumerate{
#' \item Functions which get, put, and transform data: \code{\link{nc.put.var.subset.by.axes}}, \code{\link{nc.get.var.subset.by.axes}}, \code{\link{nc.conform.data}}
#' \item Functions which deal with identifying axes, variables, and types of dimensions: \code{\link{nc.get.variable.list}}, \code{\link{nc.get.dim.axes}}, \code{\link{nc.get.dim.for.axis}}, \code{\link{nc.get.dim.bounds.var.list}}, \code{\link{nc.get.dim.names}}, \code{\link{nc.get.dim.axes.from.names}}, \code{\link{nc.get.coordinate.axes}}, \code{\link{nc.get.compress.dims}}, \code{\link{nc.is.regular.dimension}}.
#' \item Functions which deal with getting, classifying, and using time information: \code{\link{nc.get.time.series}}, \code{\link{nc.make.time.bounds}}, \code{\link{nc.get.time.multiplier}}.
#' \item Functions which make sense of projection information: \code{\link{nc.get.proj4.string}}.
#' \item Functions to ease chunked processing of data in parallel: \code{\link{get.cluster.worker.subsets}}.
#' \item Functions to ease dealing with CMIP5 data: \code{\link{get.split.filename.cmip5}}.
#' \item Utility functions: \code{\link{get.f.step.size}}.
#' }
#'
#' @name ncdf4.helpers
#' @aliases ncdf4.helpers-package
#' @docType package
#' @keywords climate ts
#' @references \url{http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions-multi.html}
NULL
