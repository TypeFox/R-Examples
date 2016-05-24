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

## these change the imputation method

setMethod("change_imputation_method", signature(data = "missing", y = "missing_variable", to = "character"), def = 
  function(y, to) {
    to <- match.arg(tolower(to) , getClass(class(y))@prototype@imputation_method)
    y@imputation_method <- to
    validObject(y, complete = TRUE)
    return(y)
  })

setMethod("change_imputation_method", signature(data = "missing", y = "missing_variable", to = "missing"), def = 
  function(y, to) {
    cat("Possible methods include:", getClass(class(y))@prototype@imputation_method, sep = "\n")
    return(invisible(NULL))
  })

setMethod("change_imputation_method", signature(data = "missing_data.frame", y = "character", to = "character"), def =
  function(data, y, to) {
    if(all(y %in% c("missing_variable", names(getClass("missing_variable")@subclasses)))) {
      mark <- sapply(colnames(data), FUN = function(x) {
        if(data@variables[[x]]@all_obs) return(FALSE)
        is(data@variables[[x]], y)
      })
      if(!any(mark)) stop(paste("no variables with missingness have class", y))
      else y <- names(mark)[mark]
    }
    y <- match.arg(y, colnames(data), several.ok = TRUE)
    if(length(to) == 1) to <- rep(to, length(y))
    else if(length(to) != length(y)) stop("'y' and 'to' must have the same length")
    for(i in 1:length(y)) {
      data@variables[[y[i]]] <- change_imputation_method(y = data@variables[[y[i]]], to = to[i])
    }
    return(data)
  })

setMethod("change_imputation_method", signature(data = "missing_data.frame", y = "numeric", to = "character"), def =
  function(data, y, to) {
    if(length(to) == 1) to <- rep(to, length(y))
    else if(length(to) != length(y)) stop("'y' and 'to' must have the same length")
    for(i in 1:length(y)) {
      data@variables[[y]] <- change_imputation_method(y = data@variables[[y]], to = to[i])
    }
    return(data)
  })

setMethod("change_imputation_method", signature(data = "missing_data.frame", y = "logical", to = "character"), def =
  function(data, y, to) {
    if(length(y) != data@DIM[2]) {
      stop("the length of 'y' must equal the number of variables in 'data'")
    }
    return(change_imputation_method(data, which(y), to))
  })

