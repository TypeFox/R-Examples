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

setMethod("change", signature(data = "missing_data.frame", y = "ANY", to = "ANY", what = "character"), def =
  function(data, y, to, what, ...) {
    if(length(what) > 1) stop("'what' must have length one")
    
    #   what <- match.arg(what, c("family", "imputation_method", "model", "size", "transformation", 
    #                             "type", "class", "link", "method"))
    if(what == "class")  what <- "type"
    if(what == "method") what <- "imputation_method"
    
    if(is.character(y) && !(what %in% c("family", "imputation_method", "link", "model", "size", "transformation", "type"))) {
      if(length(y) > 1)  stop("'y' must have length one")
      if(length(to) > 1) stop("'to' must have length one")
      if(is.logical(y) | is.numeric(y)) y <- colnames(data)[y]
      
      if(to == "unpossible") {
        mv <- data@variables[[y]]
        unpossible <- which(mv@raw_data == what)
        mv@n_unpossible <- length(unpossible)
        mv@which_unpossible <- unpossible
        
        mv@which_obs <- mv@which_obs[!(mv@which_obs %in% mv@which_unpossible)]
        mv@n_obs <- length(mv@which_obs)
        
        mv@which_miss <- mv@which_miss[!(mv@which_miss %in% mv@which_unpossible)]
        mv@n_miss <- length(mv@which_miss)
        
        data@variables[[y]] <- mv
        
        if(!length(data@weights)) {
          data@weights <- lapply(data@variables, FUN = function(y) {
            if(y@n_unpossible) {
              w <- rep(1, y@n_total)
              w[y@which_unpossible] <- 0
              return(w)
            }
            else return(NULL)
          })
        }
        else data@weights[[y]][mv@which_unpossible] <- 0
        
        return(data)
      }
      
      mv <- data@variables[[y]]
      mv@raw_data[mv@raw_data == what] <- to
      if(is.na(what) | is.na(to)) mv <- new(class(mv), raw_data = mv@raw_data, variable_name = mv@variable_name)
      data@variables[[y]] <- mv
      return(data)
    }
    
    if(what == "family") return(change_family(data = data, y = y, to = to))
    else if(what == "link") return(change_link(data = data, y = y, to = to))
    else if(what == "imputation_method") return(change_imputation_method(data = data, y = y, to = to))
    else if(what == "model") return(change_model(data = data, y = y, to = to))
    else if(what == "size") return(change_size(data = data, n = y))
    else if(what == "transformation") {
      if(missing(to)) return(change_transformation(data = data, y = y))
      else return(change_transformation(data = data, y = y, to = to, ...))
    }
    else if(what == "type") {
      if(missing(to)) return(change_type(data = data, y = y))
      else return(change_type(data = data, y = y, to = to, ...))
    }
    else stop("this should never happen")
  })

setMethod("change", signature(data = "missing_data.frame", y = "ANY", to = "numeric", what = "numeric"), def =
  function(data, y, to, what) {
    if(length(to) > 1)   stop("'to' must be a scalar")
    if(length(what) > 1) stop("'what' must be a scalar")
    if(is.logical(y) | is.numeric(y)) y <- colnames(data)[y]
    
    mv <- data@variables[[y]]
    mv@raw_data[mv@raw_data == what] <- to # NOTE: exception to "never change the raw_data slot rule"
    if(is(mv, "categorical")) {
      values <- unique(mv@raw_data)
      values <- values[!is.na(values)]
      if(length(values) == 2) mv <- new("binary", raw_data = mv@raw_data, variable_name = mv@variable_name)
    }
    else if(is.na(what) | is.na(to)) mv <- new(class(mv), raw_data = mv@raw_data, variable_name = mv@variable_name)
    else if(is(mv, "continuous")) mv@data <- mv@transformation(mv@raw_data)
    data@variables[[y]] <- mv
    return(data) ## FIXME: maybe reinitialize data?
  })

setMethod("change", signature(data = "missing_data.frame", y = "ANY", to = "logical", what = "numeric"), def =
  function(data, y, to, what) {
    change(data = data, y = y, what = what, to = as.numeric(to))
  })

setMethod("change", signature(data = "missing_data.frame", y = "ANY", to = "character", what = "numeric"), def = 
  function(data, y, to, what) {
    if(length(to) > 1)   stop("'to' must be a scalar")
    if(to != "unpossible") stop("'to' must be 'unpossible'")
    if(length(what) > 1) stop("'what' must be have length one")
    if(is.logical(y) | is.numeric(y)) y <- colnames(data)[y]
    
    mv <- data@variables[[y]]
    unpossible <- which(mv@raw_data == what)
    mv@n_unpossible <- length(unpossible)
    mv@which_unpossible <- unpossible
    
    mv@which_obs <- mv@which_obs[!(mv@which_obs %in% mv@which_unpossible)]
    mv@n_obs <- length(mv@which_obs)
    
    mv@which_miss <- mv@which_miss[!(mv@which_miss %in% mv@which_unpossible)]
    mv@n_miss <- length(mv@which_miss)
    
    data@variables[[y]] <- mv
    
    if(!length(data@weights)) {
      data@weights <- lapply(data@variables, FUN = function(y) {
        if(y@n_unpossible) {
          w <- rep(1, y@n_total)
          w[y@which_unpossible] <- 0
          return(w)
        }
        else return(NULL)
      })
    }
    else data@weights[[y]][mv@which_unpossible] <- 0
    return(data)
  })

setMethod("change", signature(data = "missing_data.frame", y = "ANY", to = "logical", what = "character"), def =
  function(data, y, to, what) {
    change(data = data, y = y, what = what, to = as.numeric(to))
  })

setMethod("change", signature(data = "mdf_list", y = "ANY", to = "ANY", what = "ANY"), def =
  function(data, y, to, what, ...) {
    out <- lapply(data, FUN = change, y = y, to = to, what = what, ...)
    class(out) <- "mdf_list"
    return(out)
  })

