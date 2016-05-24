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


setMethod("change_size", signature(data = "missing", y = "missing_variable", to = "integer"), def = 
  function(y, to) {
    n <- to
    if(n <= 0) {
      y@data <- y@data[-y@which_extra]
      y@which_extra <- integer(0)
      y@n_total <- y@n_total - y@n_extra
      y@n_extra <- NA_integer_
      return(y)
    }
    end <- y@n_total
    SEQ <- (end+1):(end+n)
    y@data <- c(y@data, rep(NA, n))
    y@which_extra <- c(y@which_extra, SEQ)
    y@n_extra <- y@n_extra + n
    y@n_total <- y@n_total + n
    return(y)
  })

setMethod("change_size", signature(data = "missing", y = "categorical", to = "integer"), def = 
  function(y, to) {
    n <- to
    if(n <= 0) {
      y@data <- y@data[-y@which_extra]
      y@which_extra <- integer(0)
      y@n_total <- y@n_total - y@n_extra
      y@n_extra <- NA_integer_
      return(y)
    }
    end <- y@n_total
    SEQ <- (end+1):(end+n)
    y@data <- c(y@data, rep(NA, n))
    y@which_extra <- c(y@which_extra, SEQ)
    y@n_extra <- y@n_extra + n
    y@n_total <- y@n_total + n
    return(y)
  })

setMethod("change_size", signature(data = "missing", y = "fixed", to = "integer"), def = 
  function(y, to) {
    n <- to
    if(n <= 0) {
      y@data <- y@data[-y@which_extra]
      y@which_extra <- integer(0)
      y@n_total <- y@n_total - y@n_extra
      y@n_extra <- NA_integer_
      return(y)
    }
    end <- y@n_total
    SEQ <- (end+1):(end+n)
    y@data <- c(y@data, rep(y@data[1], n))
    y@which_extra <- c(y@which_extra, SEQ)
    y@n_extra <- y@n_extra + n
    y@n_total <- y@n_total + n
    return(y)
  })

setMethod("change_size", signature(data = "missing_data.frame", y = "missing", to = "integer"), def = 
  function(data, to) {
    n <- to
    data@variables <- lapply(data@variables, FUN = function(x) change_size(x, n))
    data@DIM[1] <- data@variables[[1]]@n_total
    return(data)
  })
