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

## these change the transformation and inverse_transformation slots of a continuous variable

setMethod("change_transformation", signature(data = "missing", y = "missing_variable", to = "function"), def = 
  function(y, to, inverse = FALSE) {
    if(!is(y, "continuous")) stop(paste(y@variable_name, "is not a continuous variable and hence has no transformation"))
    else if(is(y, "SC_proportion")) stop(paste(y@variable_name, "is a SC_proportion and cannot change its transformation (yet)"))
    if(inverse) {
      if(identical(to, .standardize_transform)) {
        formals(to)$mean <- mean(y@raw_data, na.rm = TRUE)
        formals(to)$sd <- sd(y@raw_data, na.rm = TRUE)
      }
      else if(identical(to, .logshift)) {
        yy <- y@raw_data
        if(any(yy < 0, na.rm = TRUE)) a <- - min(yy, na.rm = TRUE)
        else a <- 0
        a <- (a + min(yy[yy > 0], na.rm = TRUE)) / 2
        formals(to)$a <- a
      }
      if("inverse" %in% names(formals(to))) formals(to)$inverse <- TRUE
      y@inverse_transformation <- to
    }
    else {
      if(identical(to, .standardize_transform)) {
        formals(to)$mean <- mean(y@raw_data, na.rm = TRUE)
        formals(to)$sd <- sd(y@raw_data, na.rm = TRUE)
      }
      else if(identical(to, .logshift)) {
        yy <- y@raw_data
        if(any(yy < 0, na.rm = TRUE)) a <- - min(yy, na.rm = TRUE)
        else a <- 0
        a <- (a + min(yy[yy > 0], na.rm = TRUE)) / 2
        formals(to)$a <- a
      }
      y@transformation <- to
      y@data <- y@transformation(y@raw_data)
    }
    return(y)
  })

setMethod("change_transformation", signature(data = "missing", y = "missing_variable", to = "missing"), def =
  function(y) {
    if(is(y, "continuous")) cat("Likely choices include:", y@known_transformations, sep = "\n")
    else cat("No transformation possible for non-continuous variables\n")
    return(invisible(NULL))
  })

setMethod("change_transformation", signature(data = "missing_data.frame", y = "character", to = "missing"), def =
  function(data, y) {
    if(all(y %in% c("continuous", names(getClass("continuous")@subclasses)))) {
      classes <- sapply(data@variables, class)
      y <- c(sapply(y, FUN = function(x) {
        names(classes[which(classes == x)])
      }))
      if(is.list(y)) stop(paste("no variables of class", names(y)[1]))
      else y <- y[1]
    }
    y <- match.arg(y, data@DIMNAMES[[2]], several.ok = TRUE)
    for(i in 1:length(y)) change_transformation(y = data@variables[[y[i]]])
    return(data)
  })

setMethod("change_transformation", signature(data = "missing_data.frame", y = "character", to = "character"), def =
  function(data, y, to) {
    if(length(to) == 1) to <- rep(to, length(y))
    else if(length(to) != length(y)) stop("'y' and 'to' must have the same length")
    if(all(y %in% c("continuous", names(getClass("continuous")@subclasses)))) {
      classes <- sapply(data@variables, class)
      y <- c(sapply(y, FUN = function(x) {
        names(classes[which(classes == x)])
      }))
      to <- rep(to[1], length(y))
    }
    y <- match.arg(y, data@DIMNAMES[[2]], several.ok = TRUE)
    trans <- lapply(to, FUN = function(x) {
      switch(x, 
             "identity" = .identity_transform,
             "standardize" = .standardize_transform,
             "squeeze" = .squeeze_transform,
             "logshift" = .logshift,
             "log" = log,
             "sqrt" = sqrt,
             "cuberoot" = .cuberoot,
             function(...) stop(paste("must replace the transformation slot for", x)))
    })
    inverse <- lapply(to, FUN = function(x) {
      switch(x, 
             "identity" = .identity_transform,
             "standardize" = .standardize_transform,
             "squeeze" = .squeeze_transform,
             "logshift" = .logshift,
             "log" = exp,
             "sqrt" = function(y, ...) y^2,
             "cuberoot" = .cuberoot,
             function(...) stop(paste("must replace the inverse_transformation slot for", x)))
    })
    for(i in 1:length(y)) {
      data@variables[[y[i]]] <- change_transformation(y = data@variables[[y[i]]], to = trans[[i]])
      data@variables[[y[i]]] <- change_transformation(y = data@variables[[y[i]]], to = inverse[[i]], inverse = TRUE)
      
      mark <- data@index[[y[i]]][1]
      data@X[,mark] <- data@variables[[y[i]]]@data
    }
    #   initialize(data)
    return(data)
  })

setMethod("change_transformation", signature(data = "missing_data.frame", y = "numeric", to = "character"), def =
  function(data, y, to) {
    return(change_transformation(data = data, y = colnames(data)[y], to = to))
  })

setMethod("change_transformation", signature(data = "missing_data.frame", y = "logical", to = "character"), def =
  function(data, y, to) {
    if(length(y) != data@DIM[2]) {
      stop("the length of 'y' must equal the number of variables in 'data'")
    }
    return(change_transformation(data = data, y = names(data@variables)[y], to = to))
  })

setMethod("change_transformation", signature(data = "missing_data.frame", y = "character", to = "function"), def =
  function(data, y, to, inverse = stop("you must specify 'inverse = FALSE' or 'inverse = TRUE'")) {
    if(all(y %in% c("continuous", names(getClass("continuous")@subclasses)))) {
      classes <- sapply(data@variables, class)
      y <- c(sapply(y, FUN = function(x) {
        names(classes[which(classes == x)])
      }))
    }
    y <- match.arg(y, data@DIMNAMES[[2]], several.ok = TRUE)
    for(i in 1:length(y)) {
      if(inverse) data@variables[[y[i]]] <- change_transformation(y = data@variables[[y[i]]], to = to, inverse = TRUE)
      else data@variables[[y[i]]] <- change_transformation(y = data@variables[[y[i]]], to = to, inverse = FALSE)
      
      mark <- data@index[[y[i]]][1]
      data@X[,mark] <- data@variables[[y[i]]]@data
    }
    return(data)
  })

setMethod("change_transformation", signature(data = "missing_data.frame", y = "numeric", to = "function"), def =
  function(data, y, to, inverse) {
    y <- names(data@variables)[y]
    return(change_transformation(data = data, y = y, to = to, inverse))
  })

setMethod("change_transformation", signature(data = "missing_data.frame", y = "logical", to = "function"), def =
  function(data, y, to, inverse) {
    if(length(y) != data@DIM[2]) {
      stop("the length of 'y' must equal the number of variables in 'data'")
    }
    return(change_transformation(data = data, y = names(data@variables)[y], to = to, inverse = inverse))
  })
