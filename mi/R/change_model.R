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

## these are convience functions that implicitly change something else by changing the model buzzword

setMethod("change_model", signature(data = "missing", y = "missing_variable", to = "character"), def = 
  function(y, to) {
    switch(to, 
           "logit" =  new("binary", variable_name = y@variable_name, raw_data = y@raw_data,
                          imputation_method = y@imputation_method, family = binomial(link = "logit")), 
           "probit" =  new("binary", variable_name = y@variable_name, raw_data = y@raw_data, 
                           imputation_method = y@imputation_method, family = binomial(link = "probit")),
           "cauchit" = new("binary", variable_name = y@variable_name, raw_data = y@raw_data, 
                           imputation_method = y@imputation_method, family = binomial(link = "cauchit")), 
           "cloglog" = new("binary", variable_name = y@variable_name, raw_data = y@raw_data, 
                           imputation_method = y@imputation_method, family = binomial(link = "cloglog")), 
           
           "qlogit" =  new("binary", variable_name = y@variable_name, raw_data = y@raw_data,
                           imputation_method = y@imputation_method, family = quasibinomial(link = "logit")), 
           "qprobit" =  new("binary", variable_name = y@variable_name, raw_data = y@raw_data, 
                            imputation_method = y@imputation_method, family = quasibinomial(link = "probit")),
           "qcauchit" = new("binary", variable_name = y@variable_name, raw_data = y@raw_data, 
                            imputation_method = y@imputation_method, family = quasibinomial(link = "cauchit")),
           "qcloglog" = new("binary", variable_name = y@variable_name, raw_data = y@raw_data, 
                            imputation_method = y@imputation_method, family = quasibinomial(link = "cloglog")),
           
           "ologit" =  new("ordered-categorical", variable_name = y@variable_name, raw_data = y@raw_data, 
                           imputation_method = y@imputation_method, family = multinomial(link = "logit")), 
           "oprobit" =  new("ordered-categorical", variable_name = y@variable_name, raw_data = y@raw_data,
                            imputation_method = y@imputation_method, family = multinomial(link = "probit")),
           "ocauchit" =  new("ordered-categorical", variable_name = y@variable_name, raw_data = y@raw_data,
                             imputation_method = y@imputation_method, family = multinomial(link = "cauchit")), 
           "ocloglog" =  new("ordered-categorical", variable_name = y@variable_name, raw_data = y@raw_data,
                             imputation_method = y@imputation_method, family = multinomial(link = "cloglog")),
           
           "mlogit" =  new("unordered-categorical", variable_name = y@variable_name, raw_data = y@raw_data,
                           imputation_method = y@imputation_method, family = multinomial(link = "logit")),
           "RNL" = new("unordered-categorical", variable_name = y@variable_name, raw_data = y@raw_data,
                       imputation_method = y@imputation_method, family = binomial(link = "logit")),
           
           "qpoisson" =  new("count", variable_name = y@variable_name, raw_data = y@raw_data, 
                             imputation_method = y@imputation_method, family = quasipoisson(link = "log")),
           "poisson" =  new("count", variable_name = y@variable_name, raw_data = y@raw_data, 
                            imputation_method = y@imputation_method, family = poisson(link = "log")),
           
           "linear" = new("continuous", variable_name = y@variable_name, raw_data = y@raw_data,
                          imputation_method = y@imputation_method, family = gaussian(link = "identity")),
           stop("model not recognized")
           )
  })

setMethod("change_model", signature(data = "missing_data.frame", y = "character", to = "character"), def =
  function(data, y, to) {
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
    check <- FALSE
    for(i in 1:length(y)) {
      categorical <- is(data@variables[[y[i]]], "categorical")
      data@variables[[y[i]]] <- change_model(y = data@variables[[y[i]]], to = to[i])
      if(categorical & !is(data@variables[[y[i]]], "categorical")) check <- TRUE
      if(!categorical & is(data@variables[[y[i]]], "categorical")) check <- TRUE
    }
    if(check) return(new(class(data), variables = data@variables))
    else      return(data)
  })

setMethod("change_model", signature(data = "missing_data.frame", y = "numeric", to = "character"), def =
  function(data, y, to) {
    if(length(to) == 1) to <- rep(to, length(y))
    else if(length(to) != length(y)) stop("'y' and 'to' must have the same length")
    for(i in 1:length(y)) {
      categorical <- is(data@variables[[y[i]]], "categorical")
      data@variables[[y[i]]] <- change_model(y = data@variables[[y[i]]], to = to[[i]])
      if(categorical & !is(data@variables[[y[i]]], "categorical")) check <- TRUE
      if(!categorical & is(data@variables[[y[i]]], "categorical")) check <- TRUE
    }
    if(check) return(new(class(data), variables = data@variables))
    else      return(data)
  })

setMethod("change_model", signature(data = "missing_data.frame", y = "logical", to = "character"), def =
  function(data, y, to) {
    if(length(y) != data@DIM[2]) {
      stop("the length of 'y' must equal the number of variables in 'data'")
    }
    return(change_model(data, which(y), to))
  })
