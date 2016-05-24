#  sysBiolAlg_fvClass.R
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
#                    definition of the class sysBiolAlg_fv                     #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_fv",
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_fv
setMethod(f = "initialize",
          signature = "sysBiolAlg_fv",
          definition = function(.Object,
                                model,
                                percentage = 100,
                                Zopt = NULL,
                                fixObjVal = TRUE,
                                tol = SYBIL_SETTINGS("TOLERANCE"),
                                lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                writeProbToFileName = NULL, ...) {

              if ( ! missing(model) ) {

                  stopifnot(is(model, "modelorg"),
                            (is.null(Zopt) || is(Zopt, "numeric")),
                            is(tol, "numeric"),
                            is(percentage, "numeric"),
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
                              paste("FV", mod_id(model), sep = "_"),
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

                  .Object <- callNextMethod(.Object,
                                            sbalg      = "fv",
                                            pType      = "lp",
                                            scaling    = scaling,
                                            fi         = 1:nCols,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = S(model),
                                            ub         = uppbnd(model),
                                            lb         = lowbnd(model),
                                            obj        = rep(0, nCols),
                                            rlb        = rep(0, nRows),
                                            rtype      = rep("E", nRows),
                                            lpdir      = lpdir,
                                            rub        = NULL,
                                            ctype      = NULL,
                                            cnames     = colNames,
                                            rnames     = rowNames,
                                            pname      = probName,
                                            algPar     = list("percentage" = percentage,
                                                              "Zopt"       = Zopt),
                                            ...)

                  # objective value
                  if ( (isTRUE(fixObjVal)) && (any(obj_coef(model) != 0)) ) {
                      if (is.null(Zopt)) {
                          optimal <- optimizeProb(model,
                                                  retOptSol = FALSE,
                                                  algorithm = "fba",
                                                  lpdir = lpdir,
                                                  scaling = scaling, ...)

                          if (optimal$ok == 0) {
                              if (lpdir == "max") {
                                  obj <- .floorValues(optimal$obj,
                                                       tol = tol)*percentage/100
                              }
                              else {
                                  obj <- .ceilValues(optimal$obj,
                                                       tol = tol)*percentage/100
                              }
                          }
                          else {
                              stop("No optimal solution!")
                          }
                      }
                      else {
                          #obj <- Zopt
                          obj <- Zopt * (percentage/100)
                      }

                      # add a row to the problem
                      #type <- ifelse(lpdir == "max", "L", "U")
                      if (lpdir == "max") {
                          type <- "L"
                          lowb <- obj
                          uppb <- SYBIL_SETTINGS("MAXIMUM")
                      }
                      else {
                          type <- "U"
                          lowb <- SYBIL_SETTINGS("MAXIMUM") * -1
                          uppb <- obj
                      }
                      oind <- which(obj_coef(model) != 0)
                      oval <- obj_coef(model)[oind]
                      addRowsToProb(lp = problem(.Object),
                                    i = met_num(model)+1,
                                    type = type, lb = lowb, ub = uppb,
                                    cind = list(oind), nzval = list(oval),
                                    rnames = "Z")
                      .Object@nr <- .Object@nr + 1L
                      .Object@alg_par[["Zopt"]] <- obj

                  }

                  if (!is.null(writeProbToFileName)) {
                      writeProb(problem(.Object),
                                fname = as.character(writeProbToFileName))
                  }

              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
