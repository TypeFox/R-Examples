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

## these change the link function used in the imputation process

setMethod("change_link", signature(data = "missing", y = "missing_variable", to = "character"), def = 
  function(y, to) {
    fam <- do.call(y@family$family, args = list(link = to))
    y@family <- fam
    validObject(y, complete = TRUE)
    return(y)
  })

setMethod("change_link", signature(data = "missing", y = "missing_variable", to = "missing"), def = 
  function(y, to) {
    cat("Likely choices include:", y@known_links, sep = "\n")
    return(invisible(NULL))
  })

setMethod("change_link", signature(data = "missing_data.frame", y = "character", to = "character"), def =
  function(data, y, to) {
    if(length(to) == 1) to <- rep(to, length(y))
    else if(length(to) != length(y)) stop("'y' and 'to' must have the same length")
    if(all(y %in% names(getClass("missing_variable")@subclasses))) {
      classes <- sapply(data@variables, class)
      y <- c(sapply(y, FUN = function(x) {
        names(classes[which(classes == x)])
      }))
      if(is.list(y)) stop(paste("no variables of class", names(y)[1]))
      else y <- y[1]
    }
    y <- match.arg(y, data@DIMNAMES[[2]], several.ok = TRUE)
    for(i in 1:length(y)) {
      data@variables[[y[i]]] <- change_link(y = data@variables[[y[i]]], to = to[i])
    }
    return(invisible(data))
  })

setMethod("change_link", signature(data = "missing_data.frame", y = "numeric", to = "character"), def =
  function(data, y, to) {
    if(length(to) == 1) to <- rep(to, length(y))
    else if(length(to) != length(y)) stop("'y' and 'to' must have the same length")
    for(i in 1:length(y)) {
      data@variables[[y]] <- change_link(y = data@variables[[y]], to = to[i])
    }
    return(invisible(data))
  })

setMethod("change_family", signature(data = "missing_data.frame", y = "logical", to = "character"), def =
  function(data, y, to) {
    if(length(y) != data@DIM[2]) {
      stop("the length of 'y' must equal the number of variables in 'data'")
    }
    return(change_link(data, which(y), to))
  })
