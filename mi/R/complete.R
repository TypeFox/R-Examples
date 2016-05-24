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


## These functions extract completed data

setMethod("complete", signature(y = "missing_variable", m = "integer"), def = 
  function(y, m, ...) {
    out <- y@data
    if(m > 0 & y@n_drawn) out[y@which_drawn] <- as.numeric(y@imputations[m,])
    return(out)
  })

setMethod("complete", signature(y = "irrelevant", m = "integer"), def = 
  function(y, m, ...) {
    return(y@raw_data)
  })

setMethod("complete", signature(y = "categorical", m = "integer"), def = 
  function(y, m, to_factor = TRUE, ...) {
    out <- y@data
    if(m > 0 & y@n_drawn) out[y@which_drawn] <- as.numeric(y@imputations[m,])
    if(to_factor) {
      out <- factor(out, ordered = is(y, "ordered-categorical"))
      levels(out) <- y@levels
    }
    return(out)
  })

setMethod("complete", signature(y = "binary", m = "integer"), def = 
  function(y, m, to_factor = TRUE, ...) {
    out <- y@data
    if(m > 0 & y@n_drawn) out[y@which_drawn] <- as.numeric(y@imputations[m,])
    if(to_factor) {
      out <- factor(out, ordered = FALSE)
      levels(out) <- y@levels
    }
    return(out)
  })

setMethod("complete", signature(y = "continuous", m = "integer"), def = 
  function(y, m, transform = TRUE, ...) {
    out <- y@data
    if(m > 0 & y@n_drawn) out[y@which_drawn] <- as.numeric(y@imputations[m,])
    if(transform) out <- y@inverse_transformation(out)
    return(out)
  })

setMethod("complete", signature(y = "nonnegative-continuous", m = "integer"), def = 
  function(y, m, transform = TRUE, ...) {
    out <- y@data
    if(m > 0 & y@n_drawn) out[y@which_drawn] <- as.numeric(y@imputations[m,])
    if(transform) {
      out <- y@inverse_transformation(out)
      out[y@raw_data == 0] <- 0
    }
    return(out)
  })

setMethod("complete", signature(y = "SC_proportion", m = "integer"), def = 
  function(y, m, transform = TRUE, ...) {
    out <- y@data
    if(m > 0 & y@n_drawn) out[y@which_drawn] <- as.numeric(y@imputations[m,])
    if(transform) out <- y@inverse_transformation(out)
    out[y@raw_data == 0] <- 0
    out[y@raw_data == 1] <- 1
    return(out)
  })

setMethod("complete", signature(y = "missing_data.frame", m = "integer"), def =
  function(y, m, to_matrix = FALSE, include_missing = TRUE) {
    if(to_matrix) out <- sapply(y@variables, complete, m = m, to_factor = FALSE, transform = FALSE)
    else          out <- as.data.frame(lapply(y@variables, complete, m = m, to_factor = TRUE, transform = TRUE))
    
    if(is(y, "allcategorical_missing_data.frame")) {
      out <- cbind(out, latents = complete(y@latents, m = m, to_factor = !to_matrix))
    }
    if(include_missing) {
      M <- is.na(y)[,!sapply(y@variables, FUN = function(y) y@all_obs), drop = FALSE]
      colnames(M) <- paste("missing", colnames(M), sep = "_")
      out <- cbind(out, M)
    }
    return(out)
  })

setMethod("complete", signature(y = "mi", m = "numeric"), def = 
  function(y, m = length(y), to_matrix = FALSE, include_missing = TRUE) {
    stopifnot(m == as.integer(m))
    m <- as.integer(m)
    l <- length(y@data)
    draws <- sum(y@total_iters)
    
    if(length(m) > 1) out <- lapply(y@data[m], complete, m = 0L, to_matrix = to_matrix, include_missing = include_missing)
    else if(m == 1)   out <- complete(y@data[[1]], m = 0L, to_matrix = to_matrix, include_missing = include_missing) # not a list
    else if(m <= l)   out <- lapply(y@data[1:m], complete, m = 0L, to_matrix = to_matrix, include_missing = include_missing)
    else { # wants more completed datasets than chains
      quotient <- m %/% l
      remainder <- m %% l
      num <- quotient + (1:l <= remainder)
      out <- vector("list", m)
      count <- 1
      for(i in seq_along(y@data)) {
        if(num[i] == 1) {
          out[[count]] <- complete(y@data[[i]], m = 0L, to_matrix = to_matrix, include_missing = include_missing)
          count <- count + 1
        }
        else { # double-dip from a chain
          SEQ <- seq(from = ceiling(draws / 2), to = draws, length.out = num[i])
          temp <- sapply(SEQ, FUN = function(j) complete(y@data[[i]], m = as.integer(j), to_matrix = to_matrix, 
                                                         include_missing = include_missing), simplify = FALSE)
          for(j in seq_along(temp)) {
            out[[count]] <- temp[[j]]
            count <- count + 1
          }
        }
      }
    }
    return(out)
  })

setMethod("complete", signature(y = "mi", m = "missing"), def = 
  function(y, to_matrix = FALSE, include_missing = TRUE) {
    return(complete(y, m = length(y), to_matrix = to_matrix, include_missing = include_missing))
  })

setMethod("complete", signature(y = "mi_list", m = "numeric"), def = 
  function(y, m = length(y[[1]]), to_matrix = FALSE, include_missing = TRUE) {
    temp <- lapply(y, FUN = complete, m = m, to_matrix = to_matrix, include_missing = include_missing)
    dfs <- temp[[1]]
    if(length(m) == 1 && m == 1 && length(temp) > 1) for(i in 2:length(temp)) {
      dfs <- rbind(dfs, temp[[i]])
    }
    else if(length(temp) > 1) for(i in 2:length(temp)) for(j in 1:length(dfs)) {
      dfs[[j]] <- rbind(dfs[[j]], temp[[i]][[j]])
    }
    return(dfs)
  })

setMethod("complete", signature(y = "mi_list", m = "missing"), def = 
  function(y, to_matrix = FALSE, include_missing = TRUE) {
    return(complete(y, m = length(y[[1]]), to_matrix = to_matrix, include_missing = include_missing))
  })
