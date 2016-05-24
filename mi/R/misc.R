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

## like sapply but for objects of mi class
mipply <- ## FIXME: should probably be a generic function instead of poor man's S4
  function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, columnwise = TRUE, to.matrix = FALSE) {
    if(is(X, "mi_list")) {
      out <- lapply(X, mipply, ..., simplify = simplify, USE.NAMES = USE.NAMES, 
                    columnwise = columnwise, to.matrix = to.matrix)
    }
    else if(is(X, "mi")) {
      X <- complete(X, to_matrix = to.matrix)
      if(columnwise) out <- sapply(X, FUN = function(x) apply(x, 2, FUN, ...), 
                                   simplify = simplify, USE.NAMES = USE.NAMES)
      else out <- sapply(X, FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES)
    }
    else if(is(X, "mdf_list")) {
      out <- lapply(X, mipply, ..., simplify = simplify, USE.NAMES = USE.NAMES, 
                    columnwise = columnwise, to.matrix = to.matrix)
    }
    else if(is(X, "missing_data.frame")) {
      if(columnwise) out <- sapply(X, FUN = function(x) apply(x, 2, FUN, ...), 
                                   simplify = simplify, USE.NAMES = USE.NAMES)
      else out <- sapply(X, FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES)
    }
    else if(is(X, "missing_variable")) {
      out <- FUN(X@data, ...)
    }
    else if(is(X, "mi_list")) {
      out <- lapply(X, FUN = mipply, ..., simplify = simplify, USE.NAMES = USE.NAMES, 
                    columnwise = columnwise, to.matrix = to.matrix)
    }
    else stop("'X' must be of class 'mi', 'missing_data.frame', 'missing_variable', or 'mi_list'")
    return(out)
  }

## create a bugs array from an mi object
mi2BUGS <-
  function(imputations, statistic = c("moments", "imputations", "parameters")) {
    if(is(imputations, "mi_list")) return(lapply(imputations, FUN = mi2BUGS, statistic = statistic))
    else if(!is(imputations, "mi")) stop("imputations must be an object of class 'mi' or 'mi_list'")

    statistic <- match.arg(statistic)
    if(statistic == "moments") {
      iterations <- sum(imputations@total_iters)
      mark <- !imputations@data[[1]]@no_missing & 
        !sapply(imputations@data[[1]]@variables, is, class2 = "irrelevant")
      means <- lapply(1:iterations, FUN = function(m) {
        matrices <- lapply(imputations@data, FUN = complete, m = m, to_matrix = TRUE, include_missing = FALSE)
        out <- sapply(matrices, colMeans)[mark,,drop = FALSE]
        return(out)
      })
      sds <- lapply(1:iterations, FUN = function(m) {
        matrices <- lapply(imputations@data, FUN = complete, m = m, to_matrix = TRUE, include_missing = FALSE)
        out <- sapply(matrices, FUN = function(x) apply(x, 2, sd))[mark,,drop = FALSE]
        return(out)
      })
      dims <- dim(means[[1]])
      arr <- array(NA_real_, c(iterations, dims[2], 2 * dims[1]), 
                   list(NULL, NULL, c(paste("mean", rownames(means[[1]]), sep = "_"), 
                                      paste("sd",   rownames(means[[1]]), sep = "_"))))
      for(i in seq_along(means)) for(j in 1:ncol(arr)) {
        arr[i,j,   1:dims[1]] <- means[[i]][,j]
        arr[i,j,-c(1:dims[1])] <- sds[[i]][,j]
      }
    }
    else if(statistic == "imputations") {
      imp_list <- lapply(imputations@data, function(x) lapply(x@variables, function(y) y@imputations))
      n.parameters <- rapply(imp_list, ncol)
      arr <- array(NA_real_, c(sum(imputations@total_iters), length(imp_list), n.parameters)) ## FIXME: names?
      for(i in seq_along(imp_list)) arr[,i,] <- unlist(imp_list[[i]])
    }
    else arr <- get_parameters(imputations)
    return(arr) # compatible with R2WinBUGS
  }

##Outputs completed data in either Stata (.dta) format or comma-separated (.csv) format
mi2stata <- 
  function(imputations, m, file, missing.ind=FALSE, ...) {
    
  if(grepl("\\.csv$", file)) type <- "csv" 
	if(grepl("\\.dta$", file)) type <- "dta" 
	
	else if(!is(imputations, "mi")) stop("imputations must be an object of class 'mi'")
	else if(!is(file, "character")) stop("filename must be specified as a character object")	
	else if(type!="dta" & type!="csv") stop("file type must be 'dta' for stata format or 'csv' for comma-separated format")

	message("Note: after loading the data into Stata, version 11 or later, type 'mi import ice' to register the data as being multiply imputed. 
	For Stata 10 and earlier, install MIM by typing 'findit mim' and include 'mim:' as a prefix for any command using the MI data.")
	unpos <- sum(sapply(imputations@data[[1]]@variables, FUN=function(x){x@n_unpossible}))
	if (unpos>0 & !missing.ind) {
		missing.ind <- TRUE
		warning("There are legitimately skipped values in the data that were not imputed.  Including variables to indicate which missing values
		were imputed.  Values which are still missing but are not indicated are legitimate skips.")
	}
	if (unpos>0 & missing.ind) {
		warning("There are legitimately skipped values in the data that were not imputed. Values which are still missing but are not indicated are legitimate skips.")
	}

	data.list <- complete(imputations, m)
	if (missing.ind) miss.indic <- data.list[[1]][,which(!is.element(colnames(data.list[[1]]), names(imputations@data[[1]]@variables)))]
	vars <- which(is.element(colnames(data.list[[1]]), names(imputations@data[[1]]@variables)))

	stata.data <- data.list[[1]][,vars]
	stata.miss <- sapply(imputations@data[[1]]@variables, FUN=function(x){
		v <- is.element(1:x@n_total, x@which_drawn)
		return(v)
	}, simplify=TRUE)
	is.na(stata.data) <- stata.miss
	if (missing.ind) stata.data <- cbind(stata.data, miss.indic)
	stata.data$mi <- 1:nrow(stata.data); stata.data$mj <- 0	

	for(i in seq_along(data.list)){
		dl <- data.list[[i]]
		if(!missing.ind) dl <- dl[,vars]
		dl$mi <- 1:nrow(dl)
		dl$mj <- i
		stata.data <- rbind(stata.data, dl)
		}
	colnames(stata.data)[which(colnames(stata.data)=="mi")] <- "_mi"
	colnames(stata.data)[which(colnames(stata.data)=="mj")] <- "_mj"

	if(type=="dta") foreign::write.dta(stata.data, file=file, version = 7L, ...)
	else if(type=="csv") write.table(stata.data, file=file, sep=",", col.names=TRUE, row.names=FALSE)
  }

## Returns the Gelman statistic
Rhats <-
  function(imputations, statistic = c("moments", "imputations", "parameters")) {
    BUGS <- mi2BUGS(imputations, statistic)
    
    make_Rhat <- function(x) {
      m <- ncol(x)
      if(m < 2) stop("need at least 2 chains to calculate an R-hat")
      iter <- nrow(x)
      xbars <- colMeans(x)
      variances <- apply(x, MARGIN = 2:3, FUN = sd)^2
      W <- colMeans(variances)
      B <- iter * apply(xbars, MARGIN = 2, FUN = var)
      R <- sqrt( (iter - 1) / iter + 1 / iter * B / W )
      return(R)
    }
    
    if(is(imputations, "mi")) return(make_Rhat(BUGS))
    else return(sapply(BUGS, FUN = make_Rhat))
  }

## tests whether a method is the one defined in my (as opposed to a user-defined method in .GlobalEnv)
is.method_in_mi <-
  function(generic, ...) {
    method <- selectMethod(generic, signature(...))
    return(environmentName(environment(method@.Data)) == "mi")
  }

## cube root transformation
.cuberoot <-
  function(y, inverse = FALSE) {
    if(inverse) y^3
    else        y^(1/3)
  }

.parse_trans <-
  function(trans) {
    if(identical(names(formals(trans)), c("y", "mean", "sd", "inverse"))) return("standardize")
    if(identical(names(formals(trans)), c("y", "a", "inverse"))) return("logshift")
    if(identical(body(trans), body(.squeeze_transform))) return("squeeze")
    if(identical(body(trans), body(.identity_transform))) return("identity")
    if(identical(body(trans), body(log))) return("log")
    if(identical(body(trans), body(sqrt))) return("sqrt")
    if(identical(body(trans), body(.cuberoot))) return("cuberoot")
    if(identical(body(trans), body(qnorm))) return("qnorm")
    return("user-defined")
  }

.prune <- function(class) {
  classes <- names(getClass(class, where = "mi")@subclasses)
  classes <- classes[!sapply(classes, isVirtualClass, where = "mi")]
  if(!isVirtualClass(class, where = "mi")) classes <- c(class, classes)
  return(classes)
}

.possible_missing_variable <-
  function(y) { ## FIXME: update this function whenever you tweak the missing_variable tree
    
    mvs <- .prune("missing_variable")
    
    maybe <- rep(TRUE, length(mvs))
    names(maybe) <- mvs
    
    if(is.factor(y)) y <- factor(y) # to drop unused levels
    vals <- unique(y)
    vals <- sort(vals[!is.na(vals)])
    
    if(length(vals) == 1) {
      maybe[] <- FALSE
      maybe["irrelevant"] <- TRUE
      maybe[.prune("fixed")] <- TRUE
      return(maybe)
    }
    else maybe[.prune("fixed")] <- FALSE
    
    if(!all(table(y) > 1)) maybe[.prune("categorical")] <- FALSE
    
    if(length(vals) == 2) { # permit binary plus children but not other kinds of categorical
      maybe[.prune("categorical")] <- FALSE
      maybe[.prune("binary")] <- TRUE
      maybe[.prune("semi-continuous")] <- FALSE
    }
    else {
      maybe[.prune("binary")] <- FALSE
    }
    
    if(!is.numeric(vals)) {
      maybe[.prune("continuous")] <- FALSE
      maybe[.prune("count")] <- FALSE
      return(maybe)
    }
    
    if(any(vals < 0)) {
      maybe[.prune("nonnegative-continuous")] <- FALSE
      maybe[.prune("positive-continuous")] <- FALSE
      maybe[.prune("count")] <- FALSE
      return(maybe)
    }
    
    if(any(vals == 0)) maybe[.prune("positive-continuous")]    <- FALSE
    else               maybe[.prune("nonnegative-continuous")] <- FALSE # unless SC_proportion
    
    if(!any(vals < 1 && vals > 0)) {
      maybe[.prune("SC_proportion")] <- FALSE
      maybe[.prune("proportion")] <- FALSE
    }
    else if(any(vals >= 1)) {
      maybe[.prune("proportion")] <- FALSE
      if(any(vals > 1)) maybe[.prune("SC_proportion")] <- FALSE
      else              maybe[.prune("SC_proportion")] <- TRUE
    }
    
    if(any(vals != as.integer(vals))) {
      maybe[.prune("count")] <- FALSE
      maybe[.prune("categorical")] <- FALSE
    }
    
    return(maybe)
  }

.cat2dummies <-
  function(y) {
    if(!is(y, "categorical")) stop("must be a categorical variable")
    if(is(y, "binary")) out <- as.matrix(as.integer(y@data == 1))
    else {
      levels <- sort(unique(y@data))
      out <- t(sapply(y@data, FUN = function(x) as.integer(x == levels)[-1]))
    }
    return(out)
  }

setMethod("fitted", signature(object = "RNL"), def = 
  function(object, ...) {
    Pr <- sapply(object, FUN = function(m) {
      eta <- m$x %*% coef(m)
      pred <- m$family$linkinv(eta)
      return(pred)
    })
    Pr <- Pr / rowSums(Pr)
    return(Pr)
  })

setMethod("fitted", signature(object = "clogit"), def = 
  function(object, ...) {
    target <- mean(as.numeric(object$y))
    lp <- object$linear.predictors
    foo <- function(par) {
      intercept <- qlogis(par)
      mean(plogis(intercept + lp)) - target
    }
    opt <- uniroot(foo, lower = 0, upper = 1)
    return(plogis(qlogis(opt$root) + lp))
  })

# Borrowed from library(MCMCpack)
.rdirichlet <-
  function(n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- rowSums(x)
    return(x / sm)
  }

