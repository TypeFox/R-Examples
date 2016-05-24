#  sysBiolAlg_roomClass.R
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


# Original code written by Marc Andre Daxer during his bachelor thesis:
# "Analysis of Gene Defects with Mixed Integer Linear Programming" (2011
# at the Heinrich-Heine-University Duesseldorf, Dpt. for Bioinformatics).
# He wrote the package sybilROOM.


#------------------------------------------------------------------------------#
#                 definition of the class sysBiolAlg_room                      #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_room",
         representation(
             wu  = "numeric",
             wl  = "numeric",
             fnc = "integer",
             fnr = "integer"
         ),
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_room
setMethod(f = "initialize",
          signature = "sysBiolAlg_room",
          definition = function(.Object,
                                model,
                                wtflux,
                                delta = 0.03,
                                epsilon = 0.001,
                                LPvariant = FALSE,
                                absMAX = SYBIL_SETTINGS("MAXIMUM"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                writeProbToFileName = NULL, ...) {

              if ( ! missing(model) ) {

                  if (missing(wtflux)) {
                      tmp <- .generateWT(model, ...)
                      wtflux <- tmp$fluxes[tmp$fldind]
                  }

                  stopifnot(is(model, "modelorg"),
                            is(wtflux, "numeric"),
                            is(delta, "numeric"),
                            is(epsilon, "numeric"),
                            is(LPvariant, "logical"),
                            is(absMAX, "numeric"))
                  
                  stopifnot(length(wtflux) == react_num(model))

                  #    the problem: minimize
                  #
                  #             |       |
                  #         S   |   0   |  = 0
                  #             |       |
                  #       -------------------------
                  #       1     |       |
                  #         1   | -vwl  |  >= wl
                  #           1 |       |
                  #       -------------------------
                  #       1     |       |
                  #         1   | -vwu  |  <= wu
                  #           1 |       |
                  #       -------------------------
                  #
                  #        lb   |   0
                  #        ub   |   1
                  # ctype: C    |   B
                  # obj:   0    |   1


                  # ---------------------------------------------
                  # problem dimensions
                  # ---------------------------------------------

                  nc <- react_num(model)
                  nr <- met_num(model)

                  nCols <- (2 * nc)
                  nRows <- (nr + 2 * nc)


                  # ---------------------------------------------
                  # constraint matrix
                  # ---------------------------------------------

                  # ROOM-Boundaries
                  
                  # if we use the formulation as linear program, delta and
                  # epsilon are set to zero (see Shlomi et al.)
                  #if (isTRUE(LPvariant)) {
                  #    wu <- wtflux
                  #    wl <- wtflux
                  #}
                  #else {
                      wu <- wtflux + delta * abs(wtflux) + epsilon
                      wl <- wtflux - delta * abs(wtflux) - epsilon
                  #}

                  vwu <- uppbnd(model) - wu
                  vwl <- lowbnd(model) - wl


                  # the initial matrix dimensions
                  LHS <- Matrix::Matrix(0, 
                                        nrow = nRows,
                                        ncol = nCols,
                                        sparse = TRUE)

                  # rows for the flux variables
                  LHS[1:nr,1:nc] <- S(model)

                  # location of the wild type strain
                  fi <- c(1:nc)

                  # constraint-matrix for ROOM
                  diag(LHS[(nr+1)    :(nr+nc),1     :nc   ]) <-        1
                  diag(LHS[(nr+nc+1) :nRows  ,1     :nc   ]) <-        1
                  diag(LHS[(nr+1)    :(nr+nc),(nc+1):nCols]) <- vwl * -1
                  diag(LHS[(nr+nc+1) :nRows  ,(nc+1):nCols]) <- vwu * -1


                  # ---------------------------------------------
                  # columns
                  # ---------------------------------------------

                  clower <- c(lowbnd(model), rep(0, nc))
                  cupper <- c(uppbnd(model), rep(1, nc))

                  if (isTRUE(LPvariant)) {
                      ctype  <- NULL
                      pt     <- "lp"
                  }
                  else {
                      ctype  <- c(rep("C", nc),  rep("B", nc))
                      pt     <- "mip"
                  }

                  # ---------------------------------------------
                  # rows
                  # ---------------------------------------------

                  #rlower <- c(rhs(model), wl, wu)
                  #rupper <- c(rhs(model), rep(absMAX, nc), wu)
                  rlower <- c(rep(0, nr), wl, wu)
                  rupper <- c(rep(0, nr), rep(absMAX, nc), wu)
                  rtype  <- c(rep("E", nr), rep("L", nc), rep("U", nc))


                  # ---------------------------------------------
                  # objective function
                  # ---------------------------------------------

                  cobj <- c(rep(0, nc), rep(1, nc))


                  # ---------------------------------------------
                  # row and column names for the problem object
                  # ---------------------------------------------

                  if (isTRUE(useNames)) {
                      if (is.null(cnames)) {
                          cn <- c(react_id(model),
                                  paste("oo", react_id(model), sep = "_")
                          )
                          colNames <- .makeLPcompatible(cn,
                                                                prefix = "x")
                      }
                      else {
                          stopifnot(is(cnames, "character"),
                                    length(cnames) == nCols)
                          colNames <- cnames
                      }

                      if (is.null(rnames)) {
                          rn <- c(met_id(model),
                                  paste("wl", react_id(model), sep = "_"),
                                  paste("wu", react_id(model), sep = "_")
                          )
                          rowNames <- .makeLPcompatible(rn,
                                                                prefix = "r")
                      }
                      else {
                          stopifnot(is(rnames, "character"),
                                    length(rnames) == nRows)
                          rowNames <- rnames
                      }

                      if (is.null(pname)) {
                          probName <- .makeLPcompatible(
                           paste("ROOM", toupper(pt), mod_id(model), sep = "_"),
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
                                            sbalg      = "room",
                                            pType      = pt,
                                            scaling    = scaling,
                                            fi         = fi,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = LHS,
                                            ub         = cupper,
                                            lb         = clower,
                                            obj        = cobj,
                                            rlb        = rlower,
                                            rub        = rupper,
                                            rtype      = rtype,
                                            ctype      = ctype,
                                            lpdir      = "min",
                                            cnames     = colNames,
                                            rnames     = rowNames,
                                            pname      = probName,
                                            algPar     = list("wtflux"  = wtflux,
                                                              "delta"   = delta,
                                                              "epsilon" = epsilon),
                                            ...)

                  .Object@wu  <- as.numeric(wu)
                  .Object@wl  <- as.numeric(wl)
                  .Object@fnr <- as.integer(nr)
                  .Object@fnc <- as.integer(nc)

                  if (!is.null(writeProbToFileName)) {
                      writeProb(problem(.Object),
                                fname = as.character(writeProbToFileName))
                  }

#                  # ---------------------------------------------
#                  # build problem object
#                  # ---------------------------------------------
#
#                  lp <- optObj(solver = solver, method = method, pType = "mip")
#                  lp <- initProb(lp, nrows = nRows, ncols = nCols)
#
#                  # ---------------------------------------------
#                  # set control parameters
#                  # ---------------------------------------------
#
#                  if (!any(is.na(solverParm))) {
#                      setSolverParm(lp, solverParm)
#                  }
#    
#
#                  loadLPprob(lp,
#                             nCols = nCols,
#                             nRows = nRows,
#                             mat   = LHS,
#                             ub    = cupper,
#                             lb    = clower,
#                             obj   = cobj,
#                             rlb   = rlower,
#                             rub   = rupper,
#                             rtype = rtype,
#                             ctype = ctype,
#                             lpdir = "min"
#                  )
#                  
#                  if (!is.null(scaling)) {
#                      scaleProb(lp, scaling)
#                  }
#
#                  .Object@problem   <- lp
#                  .Object@algorithm <- "room"
#                  .Object@nr        <- as.integer(nRows)
#                  .Object@nc        <- as.integer(nCols)
#                  .Object@fldind    <- as.integer(fi)
#                  .Object@wu        <- as.numeric(wu)
#                  .Object@wl        <- as.numeric(wl)
#                  .Object@fnr       <- as.integer(nr)
#                  .Object@fnc       <- as.integer(nc)
#                  validObject(.Object)
                  
              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                                other methods                                 #
#------------------------------------------------------------------------------#

setMethod("applyChanges", signature(object = "sysBiolAlg_room"),
    function(object, del, obj, ld,
             react    = NULL,
             lb       = NULL,
             ub       = NULL,
             obj_coef = NULL,
             fldind   = TRUE,
             lpdir    = NULL) {

        if (!isTRUE(fldind)) {
            warning("argument ", sQuote("fldind"), " is currently not used for ROOM")
        }

        fi <- fldind(object)[react]

        if (any(is.na(fi))) {
            stop("argument ", sQuote("react"), " must contain reactions only")
        }

        tmp_val <- list("fi" = react, "lb" = NULL, "ub" = NULL)

        wu    <- object@wu
        wl    <- object@wl
        lpmod <- problem(object)

        if (isTRUE(del)) {
            # store default lower and upper bounds
            tmp_val[["lb"]] <- getColsLowBnds(lpmod, fi)
            tmp_val[["ub"]] <- getColsUppBnds(lpmod, fi)
    
            # change bounds of fluxes in react
            check <- changeColsBnds(lpmod, fi, lb, ub)

            # change constraint matrix and objective coefficients
            vwu <- (ub - wu[fi]) * -1
            vwl <- (lb - wl[fi]) * -1

            ri <- react + object@fnr
            ci <- react + object@fnc
            for (i in seq(along = react)) {
                changeMatrixRow(lpmod, ri[i], c(react[i], ci[i]), c(1, vwl[i]))
                changeMatrixRow(lpmod,
                                ri[i]+object@fnc,
                                c(react[i], ci[i]),
                                c(1, vwu[i]))
            }
            
            #changeObjCoefs(lpmod, ci, rep(0, length(react)))

            tmp_val[["ri"]] <- ri
            tmp_val[["ci"]] <- ci
        }
    
        return(tmp_val)
    }
)


#------------------------------------------------------------------------------#

setMethod("resetChanges", signature(object = "sysBiolAlg_room"),
    function(object, old_val) {

        fi    <- fldind(object)
        wu    <- object@wu
        wl    <- object@wl
        lpmod <- problem(object)

        if ( (!is.null(old_val[["lb"]])) || (!is.null(old_val[["ub"]])) ) {
            check <- changeColsBnds(lpmod,
                                    fi[old_val[["fi"]]],
                                    old_val[["lb"]], old_val[["ub"]])

            vwu <- (old_val[["ub"]] - wu[fi[old_val[["react"]]]]) * -1
            vwl <- (old_val[["lb"]] - wl[fi[old_val[["react"]]]]) * -1

            for (i in seq(along = old_val[["fi"]])) {
                changeMatrixRow(lpmod,
                                old_val[["ri"]][i],
                                c(old_val[["fi"]][i], old_val[["ci"]][i]),
                                c(1, vwl[i]))
                changeMatrixRow(lpmod,
                                old_val[["ri"]][i]+object@fnc,
                                c(old_val[["fi"]][i], old_val[["ci"]][i]),
                                c(1, vwu[i]))
            }
            
            #changeObjCoefs(lpmod, old_val[["ci"]], rep(1, length(react)))
        }
    
        return(invisible(TRUE))
    }
)
