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

.guess_type <-
  function(y, favor_ordered = TRUE, favor_positive = FALSE, threshold = 5,
           variable_name = deparse(substitute(y))) {
    
    if(!is.null(dim(y))) stop(paste(variable_name, ": must be a vector"))
    if(is.factor(y)) y <- factor(y) # to drop unused levels    
    values <- unique(y)
    values <- sort(values[!is.na(values)])
    len <- length(values)
    if(len == 0) {
      warning(paste(variable_name, ": cannot infer variable type when all values are NA, guessing 'irrelevant'"))
      type <- "irrelevant"
    }
    else if(len == 1)             type <- "fixed"
    else if(grepl("^[[:punct:]]", 
                  variable_name)) type <- "irrelevant"
    else if(identical("id",
                      tolower(variable_name))) type <- "irrelevant"
    else if(len == 2) {
      if(!is.numeric(values))     type <- "binary"
      else if(all(values == 
        as.integer(values))) type <- "binary"
      else if(favor_positive) {
        if(all(values > 0))       type <- "positive-continuous"
        else if(all(values >= 0)) type <- "nonnegative-continuous"
        else                      type <- "continuous"
      }
      else                        type <- "continuous"
    }
    else if(is.ts(y)) {
      if(favor_positive) {
        if(all(values > 0))       type <- "positive-continuous"
        else if(all(values >= 0)) type <- "nonnegative-continuous"
        else                      type <- "continuous"
      }
      else                        type <- "continuous"
    }
    else if(is.ordered(y))        type <-   "ordered-categorical"
    else if(is.factor(y))         type <- "unordered-categorical"
    else if(is.character(y))      type <- "unordered-categorical"
    else if(is.numeric(y)) {
      if(all(values >= 0) && 
        all(values <= 1)) {
        
        if(any(values %in% 0:1))  type <- "SC_proportion"
        else                      type <- "proportion"
      }
      else if(len <= threshold && 
        all(values == as.integer(values)))
        type <- if(favor_ordered) "ordered-categorical" else "unordered-categorical"
      else if(favor_positive) {
        if(all(values > 0))       type <- "positive-continuous"
        else if(all(values >= 0)) type <- "nonnegative-continuous"
        else                      type <- "continuous"
      }
      else                        type <- "continuous"
    }
    else stop(paste("cannot infer variable type for", variable_name))
    
    return(type)
  }

## this constructor largely supplants typecast in previous versions of library(mi)
setMethod("missing_variable", signature(y = "ANY", type = "missing"), def = 
  function(y, favor_ordered = TRUE, favor_positive = FALSE, threshold = 5,
           variable_name = deparse(substitute(y))) {
    
    type <- .guess_type(y, favor_ordered, favor_positive, threshold, variable_name)
    return(missing_variable(y = y, type = type, variable_name = variable_name))
  })

setMethod("missing_variable", signature(y = "ANY", type = "character"), def = 
  function(y, type, variable_name = deparse(substitute(y)), ...) {
    return(new(type, raw_data = y, variable_name = variable_name, ...))
  })

.show_helper <-
  function(object) {
    type <- class(object)
    missingness <- object@n_miss
    meth <- object@imputation_method
    if(object@n_miss) {
      if(is.character(object@family)) {
        fam <- object@family
        link <- NA_character_
      }
      else {
        fam <- object@family$family
        link <- object@family$link
      }
    }
    else fam <- link <- NA_character_
    if(is(object, "continuous")) trans <- .parse_trans(object@transformation)
    else trans <- NA_character_
    
    df <- data.frame(type = type, missing = missingness, method = meth, 
                     family = fam, link = link, transformation = trans)
    rownames(df) <- object@variable_name
    if(is(object, "semi-continuous")) df <- rbind(df, .show_helper(object@indicator))
    return(df)
  }

setMethod("show", signature(object = "missing_variable"), def = 
  function(object) {
    df <- .show_helper(object)
    print(df)
  })

## setAs methods cause subtle problems with auto-coercion
# setAs(from = "unordered-categorical", to = "ordered-categorical", def =
# function(from) {
#   class(from) <- "ordered-categorical"
#   return(from)
# })
# 
# setAs(from = "ordered-categorical", to = "unordered-categorical", def =
# function(from) {
#   class(from) <- "unordered-categorical"
#   return(from)
# })
# 
# setAs(from = "binary", to = "unordered-categorical", def =
# function(from) {
#   stop("not possible or necessary to coerce from 'binary' to 'unordered-categorical'")
# })


# setAs(from = "binary", to = "ordered-categorical", def =
# function(from) {
#   stop("not possible or necessary to coerce from 'binary' to 'unordered-categorical'")
# })

# setAs(from = "nonnegative-continuous", to = "continuous", def = 
# function(from) {
#   mean <- mean(from@raw_data, na.rm = TRUE)
#   sd <- sd(from@raw_data, na.rm = TRUE)
#   from@transformation <- .standardize_transform
#   formals(from@transformation)$mean <- mean
#   formals(from@transformation)$sd <- sd
#   from@inverse_transformation <- .standardize_transform
#   formals(from@inverse_transformation)$mean <- mean
#   formals(from@inverse_transformation)$sd <- sd
#   formals(from@inverse_transformation)$inverse <- TRUE
#   from@data <- from@transformation(from@raw_data)
#   class(from) <- "continuous"
#   return(from)
# })
# 
# setAs(from = "continuous", to = "positive-continuous", def = 
# function(from) {
#   from@transformation <- log
#   from@inverse_transformation <- exp
#   from@data <- from@transformation(from@raw_data)
#   class(from) <- "positive-continuous"  
#   validObject(from)
#   return(from)
# })
# 

## maybe add more methods
## NOTE: If you change something here, look also at the change_type.R file
