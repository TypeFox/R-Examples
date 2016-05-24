#  sysBiolAlg_momaClass.R
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
#                 definition of the class sysBiolAlg_moma                      #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_moma",
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_moma
setMethod(f = "initialize",
          signature = "sysBiolAlg_moma",
          definition = function(.Object,
                                model,
                                wtflux = NULL,
                                Qmat = NULL,
                                scaleDist = NULL,
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                writeProbToFileName = NULL, ...) {

              if ( ! missing(model) ) {

                  # wild type or given flux distribution
                  if (is.null(wtflux)) {
                      tmp <- .generateWT(model, ...)
                      wtsol <- tmp$fluxes[tmp$fldind]
                  }
                  else {
                      wtsol <- wtflux
                  }

                  stopifnot(is(model, "modelorg"),
                            is(wtsol, "numeric"),
                            length(wtsol) == react_num(model),
                            is.null(scaleDist) ||
                            length(scaleDist) == react_num(model))
                  
                  #  the problem: minimize
                  #
                  #            |     
                  #         S  |  = 0
                  #            |     
                  #       -----------
                  #     lb_del
                  #     ub_del
                  #
                  #  obj sum(v_wt - v_del)^2


                  # problem dimensions
                  nCols <- react_num(model)
                  nRows <- met_num(model)

                  if (is.null(scaleDist)) {
                      sdf <- rep(1, nCols)
                  }
                  else {
                      stopifnot(is(scaleDist, "numeric"))
                      sdf <- scaleDist
                  }

                  if (is.null(Qmat)) {
                      Q <- rep(2, nCols)
                  }
                  else {
                      Q <- Qmat
                  }

                  # scaling of particular reactions in the objective function
                  Q <- Q * sdf

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
                              paste("MOMA", mod_id(model), sep = "_"),
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


                  # ---------------------------------------------
                  # build problem object
                  # ---------------------------------------------

                  .Object <- callNextMethod(.Object,
                                            sbalg      = "moma",
                                            pType      = "qp",
                                            scaling    = scaling,
                                            fi         = 1:nCols,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = S(model),
                                            ub         = uppbnd(model),
                                            lb         = lowbnd(model),
                                            #obj        = -2 * wtsol,
                                            obj        = (-2 * wtsol) * sdf,
                                            rlb        = rep(0, nRows),
                                            rtype      = rep("E", nRows),
                                            lpdir      = "min",
                                            rub        = NULL,
                                            ctype      = NULL,
                                            cnames     = colNames,
                                            rnames     = rowNames,
                                            pname      = probName,
                                            algPar     = list("wtflux" = wtsol,
                                                              "Qmat" = Q,
                                                              "scaleDist" = sdf),
                                            ...)

                  # add quadratic part of objective function
                  loadQobj(.Object@problem, Q)
                  #loadQobj(.Object@problem, rep(2, nCols))
                  #loadQobj(.Object@problem, 2 * Diagonal(nCols))

                  if (!is.null(writeProbToFileName)) {
                      writeProb(problem(.Object),
                                fname = as.character(writeProbToFileName))
                  }


#                  # make problem object
#                  lp <- optObj(solver = solver, method = method, pType = "qp")
#                  lp <- initProb(lp, nrows = nRows, ncols = nCols)
#
#                  # set parameters
#                  if (!any(is.na(solverParm))) {
#                      setSolverParm(lp, solverParm)
#                  }
#
#                  # load problem data
#                  loadLPprob(lp,
#                             nCols = nCols,
#                             nRows = nRows,
#                             mat   = S(model),
#                             ub    = uppbnd(model),
#                             lb    = lowbnd(model),
#                             obj   = -2 * wtsol,
#                             rlb   = rep(0, nRows),
#                             rub   = NULL,
#                             rtype = rep("E", nRows),
#                             lpdir = "min"
#                  )
#
#                  # load quadratic part of the objctive function
#                  loadQobj(lp, 2 * Diagonal(nCols))
#                  
#                  # scaling
#                  if (!is.null(scaling)) {
#                      scaleProb(lp, scaling)
#                  }
#
#
#                  .Object@problem   <- lp
#                  .Object@algorithm <- "moma"
#                  .Object@nr        <- as.integer(nRows)
#                  .Object@nc        <- as.integer(nCols)
#                  .Object@fldind    <- as.integer(c(1:nCols))
#                  validObject(.Object)
                  
              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
