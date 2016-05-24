#  sysBiolAlg_fbaClass.R
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


#------------------------------------------------------------------------------#
#                   definition of the class sysBiolAlg_fba                     #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_fba",
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_fba
setMethod(f = "initialize",
          signature = "sysBiolAlg_fba",
          definition = function(.Object,
                                model,
                                lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                writeProbToFileName = NULL, ...) {

              if ( ! missing(model) ) {

                  stopifnot(is(model, "modelorg"),
                            is(lpdir, "character"))
                  
                  # problem dimensions
                  nCols <- react_num(model)
                  nRows <- met_num(model)

                  # row and column names for the problem object
                  if (isTRUE(useNames)) {
                      if (is.null(cnames)) {
                          colNames <- .makeLPcompatible(react_id(model),
                                                                prefix = "x")
                      }
                      else {
                          stopifnot(is(cnames, "character"),
                                    length(cnames) == nCols)
                          colNames <- cnames
                      }

                      if (is.null(rnames)) {
                          rowNames <- .makeLPcompatible(met_id(model),
                                                                prefix = "r")
                      }
                      else {
                          stopifnot(is(rnames, "character"),
                                    length(rnames) == nRows)
                          rowNames <- rnames
                      }

                      if (is.null(pname)) {
                          probName <- .makeLPcompatible(
                              paste("FBA", mod_id(model), sep = "_"),
                              prefix = "")
                      }
                      else {
                          stopifnot(is(pname, "character"),
                                    length(pname) == 1)
                          probName <- pname
                      }
                  }
                  else {
                      colNames <- NULL
                      rowNames <- NULL
                      probName <- NULL
                  }

                  # generate problem object
                  .Object <- callNextMethod(.Object,
                                            sbalg      = "fba",
                                            pType      = "lp",
                                            scaling    = scaling,
                                            fi         = 1:nCols,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = S(model),
                                            ub         = uppbnd(model),
                                            lb         = lowbnd(model),
                                            obj        = obj_coef(model),
                                            rlb        = rep(0, nRows),
                                            rtype      = rep("E", nRows),
                                            lpdir      = lpdir,
                                            rub        = NULL,
                                            ctype      = NULL,
                                            cnames     = colNames,
                                            rnames     = rowNames,
                                            pname      = probName,
                                            ...)

                  if (!is.null(writeProbToFileName)) {
                      writeProb(problem(.Object),
                                fname = as.character(writeProbToFileName))
                  }
#
#                  # make problem object
#                  lp <- optObj(solver = solver, method = method)
#                  lp <- initProb(lp, nrows = nRows, ncols = nCols)
#
#                  # set parameters
#                  if (!any(is.na(solverParm))) {
#                      setSolverParm(lp, solverParm)
#                  }
#
#                  loadLPprob(lp,
#                             nCols = nCols,
#                             nRows = nRows,
#                             mat   = S(model),
#                             ub    = uppbnd(model),
#                             lb    = lowbnd(model),
#                             obj   = obj_coef(model),
#                             rlb   = rep(0, nRows),
#                             rub   = NULL,
#                             rtype = rep("E", nRows),
#                             lpdir = lpdir
#                  )
#                  
#                  if (!is.null(scaling)) {
#                      scaleProb(lp, scaling)
#                  }
#                  
#                  .Object@problem   <- lp
#                  .Object@algorithm <- "fba"
#                  .Object@nr        <- as.integer(nRows)
#                  .Object@nc        <- as.integer(nCols)
#                  .Object@fldind    <- as.integer(c(1:nCols))
#                  validObject(.Object)

              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
