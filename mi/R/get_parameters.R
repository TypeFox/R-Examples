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


## these extract parameters from an estimated object

setMethod("get_parameters", signature(object = "ANY"), def =
  function(object, ...) {
    return(c(coef(object)))
  })

setOldClass("polr")
setMethod("get_parameters", signature(object = "polr"), def =
  function(object, ...) {
    return(c(coef(object), object$zeta))
  })

setOldClass("multinom")
setMethod("get_parameters", signature(object = "multinom"), def =
  function(object, ...) {
    return(c(t(coef(object))))
  })

setMethod("get_parameters", signature(object = "missing_variable"), def = 
  function(object, latest = FALSE, ...) {
    if(latest) {
      if(is.logical(latest)) {
        mark <- !apply(object@parameters, 1, FUN = function(x) any(is.na(x)))
        mark <- mark[length(mark)]
      }
      else mark <- latest
      return(object@parameters[mark,])
    }
    else return(object@parameters)
  })

setMethod("get_parameters", signature(object = "missing_data.frame"), def = 
  function(object, latest = FALSE, ...) {
    mini_list <- lapply(object@variables, get_parameters, latest = latest, ...)
    out <- matrix(NA_real_, nrow(mini_list[[1]]), ncol = 0)
    for(i in seq_along(mini_list)) out <- cbind(out, mini_list[[i]])
    return(out)
  })

setMethod("get_parameters", signature(object = "mi"), def = 
  function(object, latest = FALSE, ...) {
    mini_list <- lapply(object@data, get_parameters, latest = latest, ...)
    dims <- dim(mini_list)
    out <- array(NA_real_, c(dims[1], length(mini_list), dims[2]), 
                 dimnames = list(NULL, NULL, colnames(mini_list[[1]])))
    for(i in 1:NCOL(out)) out[,i,] <- mini_list[[i]]
    return(out)
  })

setMethod("get_parameters", signature(object = "mi_list"), def = 
  function(object, latest = FALSE, ...) {
    lapply(object, get_parameters, latest = latest, ...)
  })
