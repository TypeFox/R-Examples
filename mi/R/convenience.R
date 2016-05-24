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

## Some S3 methods for convenience

as.double.missing_variable <-
  function(x, ...) {
    stop("you must write an 'as.double' method for the", class(x), "class")
  }

as.double.categorical <-
  function(x, ...) {
    x@data
  }

as.double.continuous <-
  function(x, transformed = TRUE, ...) {
    if(transformed) x@data
    else x@inverse_transformation(x@data)
  }

as.double.count <-
  function(x, ...) {
    x@data
  }

as.double.irrelevant <-
  function(x, ...) {
    as.double(x@raw_data)
  }

as.double.missing_data.frame <-
  function(x, transformed = TRUE, ...) {
    sapply(x@variables, as.double, transformed = transformed)
  }

as.data.frame.missing_data.frame <- 
  function(x, row.names = NULL, optional = FALSE, ...) {
    as.data.frame(lapply(x@variables, FUN = function(y) y@raw_data),
                  row.names = if(is.null(row.names)) rownames(x) else row.names)
  }

dim.missing_data.frame <-
  function(x) {
    x@DIM
  }

dimnames.missing_data.frame <-
  function(x) {
    x@DIMNAMES
  }

names.missing_data.frame <-
  function(x) {
    x@DIMNAMES[[2]]
  }

dim.mi <-
  function(x) {
    if(isS4(x)) x@data[[1]]@DIM else {
    	class(x) <- "list"
    	return(dim(x))
    }
  }

dimnames.mi <-
  function(x) {
    if(isS4(x)) x@data[[1]]@DIMNAMES else {
    	class(x) <- "list"
    	return(dimnames(x))
    }
  }

names.mi <-
  function(x) {
    if(isS4(x)) x@data[[1]]@DIMNAMES[[2]] else {
    	class(x) <- "list"
    	return(names(x))
    }
  }

is.na.missing_variable <-
  function(x) {
    out <- rep(FALSE, x@n_total)
    out[x@which_miss] <- TRUE
    return(out)
  }

is.na.missing_data.frame <-
  function(x) {
    sapply(x@variables, is.na)
  }

is.na.mi <-
  function(x) {
    if(isS4(x)) is.na(x@data[[1]]) else {
    	class(x) <- "list"
    	return(is.na(x))
    }
  }

length.missing_variable <-
  function(x) {
    x@n_total
  }

length.missing_data.frame <-
  function(x) {
    ncol(x)
  }

length.mi <-
  function(x) {
    if(isS4(x)) length(x@data) else {
    	class(x) <- "list"
    	return(length(x))
    }
  }

print.mdf_list <-
  function(x ,...) {
    show(x)
  }

print.mi_list <-
  function(x, ...) {
    show(x)
  }

"[.missing_data.frame" <-
  function(x, i, j, drop = if (missing(i)) TRUE else length(j) == 1) {
    if(!missing(i)) {
      cdf <- complete(x, m = 0L)
      if(!missing(j)) return(cdf[i,j,drop = drop])
      else return(cdf[i,,drop = drop])
    }
    else if(length(j) > 1) return(new(class(x), variables = x@variables[j]))
    else if(is.numeric(j) && j < 0) return(new(class(x), variables = x@variables[j]))
    else return(x@variables[[j]])
  }

"[<-.missing_data.frame" <-
  function (x, i, j, value) {
    if(!missing(i)) {
      if(!missing(j)) x@variables[[j]]@raw_data[i,] <- value
      else stop("a variable (column) must be specified when replacing")
    }
    else if(is.null(value)) x@variables[j] <- value
    else if(is(value, "missing_variable")) x@variables[[j]] <- value
    else stop("replacement must be 'NULL' or a 'missing_variable'")
    return(new(class(x), variables = x@variables))
  }

"[[.missing_data.frame" <-
  function(x, ..., exact = TRUE) {
    return(x[,...])
  }

"[[<-.missing_data.frame" <-
  function (x, i, j, value) {
    if(missing(j)) x[,i] <- value
    else x[i,j] <- value
    return(x)
  }

"$.missing_data.frame" <-
  function(x, name) {
    return(x[,name])
  }


"$<-.missing_data.frame" <-
  function(x, name, value) { # this never gets dispatched for some reason
    x[,name] <- value
    return(x)
  }
