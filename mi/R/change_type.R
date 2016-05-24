# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## these coerce a missing_variable to a different type of missing_variable

setMethod("change_type", signature(data = "missing", y = "missing_variable", to = "character"), def = 
  function(y, to, ...) {
    to <- match.arg(to, names(getClass("missing_variable")@subclasses))
    if(to %in% c("ordered-categorical", "binary")) raw <- as.ordered(y@raw_data)
    else if(to == "unordered-categorical")         raw <- factor(y@raw_data, ordered = FALSE)
    else raw <- as.numeric(y@raw_data)
    
    vals <- unique(raw)
    vals <- vals[!is.na(vals)]
    if(length(vals) <= 1) {
      warning(paste(y@variable_name, ": cannot change type because only one unique value"))
      return(y)
    }
    else return(new(to, variable_name = y@variable_name, raw_data = raw, imputation_method = y@imputation_method, ...))
  })

setMethod("change_type", signature(data = "missing", y = "missing_variable", to = "missing"), def = 
  function(y, to) {
    classes <- .possible_missing_variable(y@raw_data)
    classes <- names(classes[classes])
    cat("Likely choices include:", classes, sep = "\n")
    return(invisible(NULL))
  })

setMethod("change_type", signature(data = "missing_data.frame", y = "character", to = "missing"), def =
  function(data, y, to) {
    if(all(y %in% names(getClass("missing_variable")@subclasses))) {
      classes <- sapply(data@variables, class)
      y <- c(sapply(y, FUN = function(x) {
        names(classes[which(classes == x)])
      }))
      if(is.list(y)) stop(paste("no variables of class", names(y)[1]))
      else y <- y[1]
    }
    y <- match.arg(y, data@DIMNAMES[[2]], several.ok = TRUE)
    for(i in 1:length(y)) change_type(y = data@variables[[y[i]]])
    return(data)
  })

setMethod("change_type", signature(data = "missing_data.frame", y = "character", to = "character"), def =
  function(data, y, to, ...) {
    if(length(to) == 1) to <- rep(to, length(y))
    else if(length(to) != length(y)) stop("'y' and 'to' must have the same length")
    if(all(y %in% names(getClass("missing_variable")@subclasses))) {
      classes <- sapply(data@variables, class)
      y <- c(sapply(y, FUN = function(x) {
        names(classes[which(classes == x)])
      }))
      if(is.list(y)) stop(paste("no variables of class", names(y)[1]))
      to <- rep(to[1], length(y))
    }
    y <- match.arg(y, data@DIMNAMES[[2]], several.ok = TRUE)
    for(i in 1:length(y)) {
      data@variables[[y[i]]] <- change_type(y = data@variables[[y[i]]], to = to[i], ...)
      data@variables[[y[i]]]@variable_name = y[i]
    }
    return(new(class(data), variables = data@variables))
  })

setMethod("change_type", signature(data = "missing_data.frame", y = "numeric", to = "character"), def =
  function(data, y, to, ...) {
    if(length(to) == 1) to <- rep(to, length(y))
    else if(length(to) != length(y)) stop("'y' and 'to' must have the same length")
    for(i in 1:length(y)) {
      data@variables[[y]] <- change_type(y = data@variables[[y]], to = to[[i]], ...)
    }
    return(new(class(data), variables = data@variables))
  })

setMethod("change_type", signature(data = "missing_data.frame", y = "logical", to = "character"), def =
  function(data, y, to, ...) {
    if(length(y) != data@DIM[2]) {
      stop("the length of 'y' must equal the number of variables in 'data'")
    }
    return(change_type(data, which(y), to, ...))
  })

