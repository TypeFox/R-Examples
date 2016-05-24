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

setMethod("missing_data.frame", signature(y = "data.frame"), def = 
  function(y, subclass = NA_character_, by = NULL, types = NULL, 
           favor_ordered = TRUE, favor_positive = FALSE, threshold = 5, ...) {
    if(!is.na(subclass) && subclass == "allcategorical") threshold <- Inf
    if(!is.null(by)) {
      mdfs <- by(y, lapply(by, FUN = function(b) y[,b]), FUN = function(d) {
        missing_data.frame(d, favor_ordered = favor_ordered, threshold = threshold,
                           favor_positive = favor_positive, subclass = subclass, types = types, ...)
      })
      class(mdfs) <- "mdf_list"
      return(mdfs)
    }
    variables <- vector("list", length = ncol(y))
    if(is.null(types)) for(i in seq_along(variables)) {
      variables[[i]] <- missing_variable(y[,i], favor_ordered = favor_ordered, favor_positive = favor_positive, 
                                         variable_name = colnames(y)[i], threshold = threshold)
    }
    else for(i in seq_along(variables)) {
      variables[[i]] <- new(types[i], variable_name = colnames(y)[i], raw_data = y[,i])
    }
    
    if(is.na(subclass)) new("missing_data.frame", variables = variables, DIMNAMES = dimnames(y), ...)
    else new(paste(subclass, "missing_data.frame", sep = "_"), variables = variables, DIMNAMES = dimnames(y), ...)
  }
          )

setMethod("missing_data.frame", signature(y = "matrix"), def = 
  function(y, ...) {
    return(missing_data.frame(y = as.data.frame(y), ...))
  }
          )

setMethod("missing_data.frame", signature(y = "list"), def = 
  function(y, ...) {
    if(!all(sapply(y, is, class2 = "missing_variable"))) {
      stop("all list elements must inherit from the 'missing_variable' class")
    }
    return(new("missing_data.frame", variables = y))
  }
          )

setAs(from = "data.frame", to = "missing_data.frame", def =
  function(from) {
    missing_data.frame(from)
  }
      )

setAs(from = "matrix", to = "missing_data.frame", def = 
  function(from) {
    missing_data.frame(as.data.frame(from))
  }
      )

setAs(from = "missing_data.frame", to = "data.frame", def =
  function(from) {
    return(complete(from, m = 0L))
  }
      )

setAs(from = "missing_data.frame", to = "matrix", def =
  function(from) {
    return(complete(from, m = 0L, to_matrix = TRUE))
  }
      )

## FIXME: Probably need to add a boatload of methods to mimic the behavior of data.frames

.default_model <-
  function(y, data) {
    if(is(data, "allcategorical_missing_data.frame")) return("Gibbs")
    if(y@all_obs) {
      if(is(y, "semi-continuous")) return(rep(NA_character_, 2))
      else return(NA_character_)
    }
    if(is(y, "irrelevant")) return(NA_character_)
    if(y@imputation_method == "mcar") return(NA_character_)
    if(!is.method_in_mi("fit_model", y = class(y), data = class(data))) {
      if(is(y, "semi-continuous")) return(rep("user-defined", 2))
      else return("user-defined")
    }
    fam <- y@family$family
    link <- y@family$link
    if(is(y, "count")) {
      if(fam == "quasipoisson" && link == "log") return("qpoisson")
      else if(fam == "poisson" && link == "log") return("poisson")
      else return("****")
    }
    else if(is(y, "binary")) {
      if(is(y, "grouped-binary")) return("clogit")
      if(fam == "quasibinomial")  return(paste("q", link, sep = ""))
      else if(fam == "binomial")  return(link)
      else return("****")
    }
    else if(is(y, "interval")) return("survreg")
    else if(is(y, "ordered-categorical"))   return(paste("o", link, sep = ""))
    else if(is(y, "unordered-categorical")) {
      if(fam == "binomial") out <- "RN"
      else out <- "m"
      return(paste(out, link, sep = ""))
    }
    else if(is(y, "proportion")) return(if(fam == "gaussian") "linear" else "betareg")
    else if(is(y, "SC_proportion")) {
      out <- .default_model(y@indicator, data)
      return(c("betareg", out))
    }
    else if(is(y, "semi-continuous")) {
      out <- .default_model(y@indicator, data)
      if(fam == "gaussian") { 
        if(link == "identity") return(c("linear", out))
        else if(link == "log") return(c("loglinear", out))
        else if(link == "inverse") return(c("inverselinear", out))
        else return(c("****", out))
      }
      else if(fam == "Gamma") return(c("****", out))
      else if(fam == "inverse.gaussian") return(c("****", out))
      else if(fam == "quasi") return(c("quasi", out))
      else return(c("****", out))
    }
    else if(is(y, "continuous")) {
      if(fam == "gaussian") { 
        if(link == "identity") return("linear")
        else if(link == "log") return("loglinear")
        else if(link == "inverse") return("inverselinear")
        else return("****")
      }
      else if(fam == "Gamma") return("****")
      else if(fam == "inverse.gaussian") return("****")
      else if(fam == "quasi") return("quasi")
      else return("****")
    }
    else return("user-defined")
  }

setMethod("show", "missing_data.frame", def = 
  function(object) {
    k <- object@DIM[2]
    df <- .show_helper(object@variables[[1]])
    for(i in 2:k) {
      df <- rbind(df, .show_helper(object@variables[[i]]))
    }
    df1 <- cbind(df[,1:3], model = unlist(sapply(object@variables, FUN = function(y) .default_model(y, object))))
    if(is(object, "experiment_missing_data.frame")) df1$concept[names(object@concept)] <- object@concept
    df2 <- df[,-c(1:3)]
    cat("Object of class", class(object), "with", nrow(object), "observations on", ncol(object), "variables\n")
    if(length(object@patterns)) {
      npatterns <- nlevels(object@patterns)
      cat("\nThere are", npatterns, "missing data patterns\n")
      #       print(table(as.integer(object@patterns)))
      #       mat <- as.matrix(levels(object@patterns))
      #       colnames(mat) <- "missing"
      #       print(mat)
      cat("\nAppend '@patterns' to this", class(object), "to access the corresponding pattern for every observation or perhaps use table()\n\n")
    }
    print(df1)
    cat("\n")
    print(df2)
    if(any(df1$model == "****", na.rm = TRUE)) {
      cat("\n**** The model lacks a widely-recognized name but is determined by the chosen type, family, and link.\n")
    }
    return(invisible(NULL))
  })

setMethod("show", "mdf_list", def = 
  function(object) {
    for(i in seq_along(object)) {
      cat("\n", names(object)[i], "\n")
      show(object[[i]])
    }
    return(invisible(NULL))
  })

setMethod("summary", "missing_data.frame", def = 
   function(object) {
     summary(as.data.frame(object))
   })
