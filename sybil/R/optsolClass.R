#  optsolClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


# optsolClass


#------------------------------------------------------------------------------#
#                            class definitions                                 #
#------------------------------------------------------------------------------#

setClass("optsol",
    representation(
        mod_id       = "character",        # model id of the original model
        mod_key      = "character",        # model key of the original model
        solver       = "character",        # the used lp solver
        method       = "character",        # the used method
        algorithm    = "character",        # the used algorithm
        num_of_prob  = "integer",          # number of problems to solve
        lp_num_cols  = "integer",          # number of reactions
        lp_num_rows  = "integer",          # number of metabolites
        lp_obj       = "numeric",          # solution of the objective function
        lp_ok        = "integer",          # exit status of the lp solver
        lp_stat      = "integer",          # solution status
        lp_dir       = "factor",           # direction of optimization(s)
        obj_coef     = "numeric",          # objective coefficients in model
        obj_func     = "character",        # objective function (printObjFunc)
        fldind       = "integer",          # indices of fluxes
        fluxdist     = "fluxDistribution", # the flux distribution
        alg_par      = "list"              # parameters to the algorithm
    ),
    contains = "VIRTUAL",
    #validity = .validoptsol
)


#------------------------------------------------------------------------------#
#                              user constructors                               #
#------------------------------------------------------------------------------#

optsol <- function(solver) {
    if (missing(solver)) {
        stop("Creating an object of class optsol needs a valid solver!")
    }
    solver <- as.character(solver)
    new("optsol", solver = solver)
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# mod_id
setMethod("mod_id", signature(object = "optsol"),
          function(object) {
              return(object@mod_id)
          }
)

setReplaceMethod("mod_id", signature = (object = "optsol"),
                 function(object, value) {
                     object@mod_id <- value
                     return(object)
                 }
)


# mod_key
setMethod("mod_key", signature(object = "optsol"),
          function(object) {
              return(object@mod_key)
          }
)

setReplaceMethod("mod_key", signature = (object = "optsol"),
                 function(object, value) {
                     object@mod_key <- value
                     return(object)
                 }
)


# solver
setMethod("solver", signature(object = "optsol"),
          function(object) {
              return(object@solver)
          }
)

setReplaceMethod("solver", signature = (object = "optsol"),
                 function(object, value) {
                     object@solver <- value
                     return(object)
                 }
)


# method
setMethod("method", signature(object = "optsol"),
          function(object) {
              return(object@method)
          }
)

setReplaceMethod("method", signature = (object = "optsol"),
                 function(object, value) {
                     object@method <- value
                     return(object)
                 }
)


# method
setMethod("algorithm", signature(object = "optsol"),
          function(object) {
              return(object@algorithm)
          }
)

setReplaceMethod("algorithm", signature = (object = "optsol"),
                 function(object, value) {
                     object@algorithm <- value
                     return(object)
                 }
)


# num_of_prob
setMethod("num_of_prob", signature(object = "optsol"),
          function(object) {
              return(object@num_of_prob)
          }
)

setReplaceMethod("num_of_prob", signature = (object = "optsol"),
                 function(object, value) {
                     object@num_of_prob <- value
                     return(object)
                 }
)


# lp_num_cols
setMethod("lp_num_cols", signature(object = "optsol"),
          function(object) {
              return(object@lp_num_cols)
          }
)

setReplaceMethod("lp_num_cols", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_num_cols <- value
                     return(object)
                 }
)


# lp_num_rows
setMethod("lp_num_rows", signature(object = "optsol"),
          function(object) {
              return(object@lp_num_rows)
          }
)

setReplaceMethod("lp_num_rows", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_num_rows <- value
                     return(object)
                 }
)


# lp_dir
setMethod("lp_dir", signature(object = "optsol"),
          function(object) {
              return(object@lp_dir)
          }
)

setReplaceMethod("lp_dir", signature(object = "optsol", value = "factor"),
                 function(object, value) {
                     object@lp_dir <- value
                     return(object)
                 }
)

setReplaceMethod("lp_dir", signature(object = "optsol", value = "character"),
                 function(object, value) {
                     if (all(value == "min" | value == "max")) {
                         object@lp_dir <- factor(value)
                     }
                     else {
                         warning("only 'min' (1) and 'max' (-1) are allowed")
                     }
                     return(object)
                 }
)

setReplaceMethod("lp_dir", signature(object = "optsol", value = "numeric"),
                 function(object, value) {
                     if (all(value == 1 | value == -1)) {
                         fc <- value
                         fc[fc ==  1] <- "min"
                         fc[fc == -1] <- "max"
                         object@lp_dir <- factor(fc)
                     }
                     else {
                         warning("only '1' (min) and '-1' (max) are allowed")
                     }
                     return(object)
                 }
)


# lp_obj
setMethod("lp_obj", signature(object = "optsol"),
          function(object) {
              return(object@lp_obj)
          }
)

setReplaceMethod("lp_obj", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_obj <- value
                     return(object)
                 }
)


# lp_ok
setMethod("lp_ok", signature(object = "optsol"),
          function(object) {
              return(object@lp_ok)
          }
)

setReplaceMethod("lp_ok", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_ok <- value
                     return(object)
                 }
)


# lp_stat
setMethod("lp_stat", signature(object = "optsol"),
          function(object) {
              return(object@lp_stat)
          }
)

setReplaceMethod("lp_stat", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_stat <- value
                     return(object)
                 }
)


# objective coefficient
setMethod("obj_coef", signature(object = "optsol"),
          function(object) {
              return(object@obj_coef)
          }
)

setReplaceMethod("obj_coef", signature(object = "optsol"),
          function(object, value) {
              object@obj_coef <- value
              return(object)
          }
)


# objective function
setMethod("obj_func", signature(object = "optsol"),
          function(object) {
              return(object@obj_func)
          }
)

setReplaceMethod("obj_func", signature(object = "optsol"),
          function(object, value) {
              object@obj_func <- value
              return(object)
          }
)


# fldind
setMethod("fldind", signature(object = "optsol"),
          function(object) {
              return(object@fldind)
          }
)

setReplaceMethod("fldind", signature = (object = "optsol"),
                 function(object, value) {
                     object@fldind <- value
                     return(object)
                 }
)


# fluxdist
setMethod("fluxdist", signature(object = "optsol"),
          function(object) {
              return(object@fluxdist)
          }
)

setReplaceMethod("fluxdist", signature = (object = "optsol"),
                 function(object, value) {
                     object@fluxdist <- value
                     return(object)
                 }
)


# fluxes
setMethod("fluxes", signature(object = "optsol"),
          function(object) {
              return(fluxes(object@fluxdist))
          }
)

setReplaceMethod("fluxes", signature = (object = "optsol"),
                 function(object, value) {
                     fluxes(object@fluxdist) <- value
                     return(object)
                 }
)


# alg_par
setMethod("alg_par", signature(object = "optsol"),
          function(object) {
              return(object@alg_par)
          }
)

setReplaceMethod("alg_par", signature = (object = "optsol"),
                 function(object, value) {
                     object@alg_par <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

# mod_obj
setMethod("mod_obj", signature(object = "optsol"),
          function(object) {
              if ( (is.na(fldind(object)[1L])) ||
                   (is.na(fluxes(object)[1L, 1L])) ) {
                  val <- lp_obj(object)
              }
              else {
                  val <- crossprod(obj_coef(object),
                                   fluxes(object)[fldind(object), ])[1L, ]
              }
              return(val)
          }
)


# number of fluxes
setMethod("nfluxes", signature(object = "optsol"),
          function(object) {
              return(num_of_fluxes(object@fluxdist))
          }
)


# check solution status
setMethod("checkStat", signature(opt = "optsol"),
          function(opt) {
              ng <- checkSolStat(opt@lp_stat, opt@solver)
              return(ng)
          }
)


# get part of the flux distribution
setMethod("getFluxDist", signature(lp = "optsol"),
          function(lp, react = NULL, opt = NULL, drop = TRUE) {
          
              if (num_of_fluxes(fluxdist(lp)) == 1) {
                  return(fluxes(lp))
              }
              
              # cr: row indices (reactions)
              # nr: row names (reaction id's, or row indices)
              # cc: column indices (optimizations)
              # nc: number of optimization
                  
              if (is.null(react)) {
                  cr <- fldind(lp)
                  nr <- cr
              }
              else {
                  if (is(react, "reactId")) {
                      stopifnot(identical(mod_id(lp), mod_id(react)))
                      
                      cr <- fldind(lp)[react_pos(react)]
                      nr <- react_id(react)
                  }
                  else if (is(react, "numeric")) {
                      if (max(react) > nvar(fluxdist(lp))) {
                          stop("react must be in [1,", nvar(fluxdist(lp)),"]")
                      }
                      else {
                          cr <- react
                          nr <- react
                      }
                  }
                  else {
                      stop("react must be numeric or of class reactId")
                  }
              }

              if (is.null(opt)) {
                  cc <- 1:ncol(fluxes(lp))
                  nc <- cc
              }
              else {
                  if (is(opt, "numeric")) {
                      if (max(opt) > ncol(fluxes(lp))) {
                          stop("opt must be in [1,", ncol(fluxes(lp)),"]")
                      }
                      else {
                          cc <- opt
                          nc <- cc
                      }
                  }
                  else {
                      stop("opt must be numeric")
                  }
              }

              
              fl <- fluxes(lp)[cr, cc, drop = drop]
              if (is.null(dim(fl))) {
                  if (length(nr) == length(fl)) {
                      names(fl) <- nr
                  }
              }
              else {
                  rownames(fl) <- nr
                  colnames(fl) <- nc
                  #colnames(fl) <- paste("[", 1:num_of_prob(lp), "]", sep = "")
              }

              return(fl)
          }
)


# consider using sprintf here
setMethod("show", signature(object = "optsol"),
    function(object) {
        cat(sprintf("%-42s%s\n", "solver:", solver(object)))
        cat(sprintf("%-42s%s\n", "method:", method(object)))
        cat(sprintf("%-42s%s\n", "algorithm:", algorithm(object)))
        cat(sprintf("%-42s%s\n", "number of variables:", lp_num_cols(object)))
        cat(sprintf("%-42s%s\n", "number of constraints:", lp_num_rows(object)))
        if (num_of_prob(object) == 1) {
            cat(sprintf("%-42s%s\n", "return value of solver:", getMeanReturn(lp_ok(object), solver(object))))
            cat(sprintf("%-42s%s\n", "solution status:", getMeanStatus(lp_stat(object), solver(object))))
            cat(sprintf("%s%-14s%f\n", "value of objective function ", paste("(", algorithm(object), "):", sep = ""), lp_obj(object)))
            cat(sprintf("%-42s%f\n", "value of objective function (model):", mod_obj(object)))
            if ( (algorithm(object) == "moma") && (!is.null(alg_par(object)[["wtflux"]])) ) {
                cat(sprintf("%-42s%f\n", "euclidean distance between wt and v:",
                            sqrt(crossprod(alg_par(object)[["wtflux"]] - getFluxDist(object)))))
            }
        }
        else {
            cat(sprintf("%-42s%s\n", "number of problems to solve:", num_of_prob(object)))
            ok <- sum(lp_ok(object) == 0, na.rm = TRUE)
            cat(sprintf("%-42s%s\n", "number of successful solution processes:", ok))
        }
    }
)


# length of an object of class optsol
setMethod("length", signature = signature(x = "optsol"),
          function(x) {
              return(num_of_prob(x))
          }
)


# draw a histogramm (package lattice)
setMethod("plot", signature(x = "optsol", y = "missing"),
          function(x, y,
                   col = "grey",
                   xlab = "value of objective function", ...) {

              histogram(x = mod_obj(x), col = col, xlab = xlab, ...)
              
          }
)


# checkOptSol
setMethod("checkOptSol", signature(opt = "optsol"),
          function(opt, onlywarn = FALSE) {

    lp_check <- FALSE

    if (isTRUE(onlywarn)) {
        if (sum(lp_ok(opt) != 0) != 0) {
            msg <- paste("some optimizations did not end successful",
                         "check results with checkOptSol()", sep = "\n")
            warning(msg, call. = FALSE)
        }
        else {
            lp_check <- TRUE
        }
    }
    else {
		lp_check <- checksol()

		if (solver(opt) == "cplexAPI") {
			TEMP_ENV <- cplexAPI::openEnvCPLEX()
		}
		else {
			TEMP_ENV <- NULL
		}

		num_of_prob(lp_check) <- length(opt)
		
		ec <- table(lp_ok(opt))
		sc <- table(lp_stat(opt))

		exit_code(lp_check)    <- as.integer(rownames(ec))
		exit_num(lp_check)     <- as.integer(ec)
		exit_meaning(lp_check) <- mapply(getMeanReturn, exit_code(lp_check),
									   MoreArgs = list(solver = solver(opt)))

		status_code(lp_check)    <- as.integer(rownames(sc))
		status_num(lp_check)     <- as.integer(sc)
		status_meaning(lp_check) <- mapply(getMeanStatus, status_code(lp_check),
										MoreArgs = list(solver = solver(opt),
										env = TEMP_ENV))

		if (!is.null(TEMP_ENV)) {
			cplexAPI::closeEnvCPLEX(TEMP_ENV)
		}
    }

    return(lp_check)
              
}
)
