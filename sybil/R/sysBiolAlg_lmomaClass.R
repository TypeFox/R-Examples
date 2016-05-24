#  sysBiolAlg_lmomaClass.R
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
#                 definition of the class sysBiolAlg_lmoma                     #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_lmoma",
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_lmoma
setMethod(f = "initialize",
          signature = "sysBiolAlg_lmoma",
          definition = function(.Object,
                                model,
                                wtflux,
                                COBRAflag = FALSE,
                                wtobj = NULL,
                                wtobjLB = TRUE,
                                obj_coefD = NULL,
                                absMAX = SYBIL_SETTINGS("MAXIMUM"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                writeProbToFileName = NULL, ...) {

              if ( ! missing(model) ) {

                  if (missing(wtflux)) {
                      tmp    <- .generateWT(model, ...)
                      wtflux <- tmp$fluxes[tmp$fldind]
                      wtobj  <- tmp$obj
                  }
                  
                  stopifnot(is(model, "modelorg"),
                            is(COBRAflag, "logical"),
                            is(wtobjLB, "logical"),
                            is(wtflux, "numeric"),
                            is(absMAX, "numeric"))
                  
                  stopifnot(length(wtflux) == react_num(model))

                  #  the problem: minimize
                  #
                  #            |      |      |                ]
                  #        Swt |  0   |  0   |  = 0           ] left out if not COBRA
                  #            |      |      |                ]
                  #       -------------------------
                  #            |      |      |
                  #         0  | Sdel |  0   |  = 0
                  #            |      |      |
                  #       -------------------------
                  #            |      |      |
                  #         v1 - v2   |delta-| >= 0
                  #         v2 - v1   |delta+| >= 0
                  #            |      |      |
                  #       -------------------------
                  #       c_wt |  0   |  0   | >= c^T * v_wt  ] left out if not COBRA
                  #            |      |      |
                  #  lb   v_wt |del_lb|  0   |
                  #  ub   v_wt |del_ub|10000 |
                  #            |      |      |
                  #            |      |      |
                  #  obj    0  |  0   |  1   |


                  # ---------------------------------------------
                  # problem dimensions
                  # ---------------------------------------------

                  nc     <- react_num(model)
                  nr     <- met_num(model)

                  nCols  <- 4*nc
                  nRows  <- ifelse(isTRUE(COBRAflag),
                                   2*nr + 2*nc + 1,
                                   nr + 2*nc)

                  if (is.null(obj_coefD)) {
                      deltaobj <- rep(1, 2*nc)
                  }
                  else {
                      stopifnot(length(obj_coefD) == 2*nc)
                      deltaobj <- obj_coefD
                  }
                  

                  # ---------------------------------------------
                  # constraint matrix
                  # ---------------------------------------------

                  # the initial matrix dimensions
                  LHS <- Matrix::Matrix(0, 
                                        nrow = nr + 2*nc,
                                        ncol = nCols,
                                        sparse = TRUE)

                  # rows for the wild type strain
                  LHS[1:nr,(nc+1):(2*nc)] <- S(model)

                  # location of the wild type strain
                  fi <- c((nc+1):(2*nc))

                  # rows for the delta match matrix
                  diag(LHS[(nr+1)   :(nr+2*nc),1       :(2*nc)]) <- -1
                  diag(LHS[(nr+1)   :(nr+2*nc),(2*nc+1):(4*nc)]) <-  1
                  diag(LHS[(nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)]) <-  1
                  diag(LHS[(nr+nc+1):(nr+2*nc),1       :nc    ]) <-  1


                  # contraint matrix for linearMOMA, COBRA version
                  if (isTRUE(COBRAflag)) {
          
                      # rows for the wild type strain
                      LHSwt <- Matrix::Matrix(0,
                                              nrow = nr,
                                              ncol = nCols,
                                              sparse = TRUE)
                      LHSwt[1:nr,1:nc] <- S(model)
          
                      # fix the value of the objective function
                      crow <- Matrix::Matrix(c(obj_coef(model), rep(0, 3*nc)),
                                             nrow = 1,
                                             ncol = nCols,
                                             sparse = TRUE)
          
                      # the final contraint matrix
                      LHS <- rBind(LHSwt, LHS, crow)
          
                      subalg <- "lmoma_cobra"
                  }
                  else {
                      subalg <- "lmoma"
                  }


                  # ---------------------------------------------
                  # lower and upper bounds
                  # ---------------------------------------------

                  if (isTRUE(COBRAflag)) {
                      # Here we calculate wild type and deletion strain
                      # simultaineously, so we need upper and lower bounds
                      # for both, the wild type and the deletion strain.
                      # All the delta's are positive.
                      lower <- c(lowbnd(model),
                                 lowbnd(model),
                                 rep(0, 2*nc))
                      upper <- c(uppbnd(model),
                                 uppbnd(model),
                                 rep(absMAX, 2*nc))
          
          
                      rlower <- c(rep(0, nRows-1), wtobj)
                      rupper <- c(rep(0, 2*nr), rep(absMAX, 2*nc), wtobj)
                      #rupper <- rlower
                      #rupper <- c(rep(0, 2*nr), rep(0, 2*nc), 0)
                  }
                  else {
                      # Here, we keep the wild type flux distribution fixed.
                      lower  <- c(wtflux, lowbnd(model), rep(0, 2*nc))
                      upper  <- c(wtflux, uppbnd(model), rep(absMAX, 2*nc))
                      rlower <- c(rep(0, nRows))
                      rupper <- c(rep(0, nr), rep(absMAX, 2*nc))
                  }

                  # ---------------------------------------------
                  # constraint type
                  # ---------------------------------------------

                  if (isTRUE(COBRAflag)) {
                      rtype <- c(rep("E", 2*nr), rep("L", 2*nc))
                      if (isTRUE(wtobjLB)) {
                          rtype <- append(rtype, "L")
                      }
                      else {
                          rtype <- append(rtype, "U")
                      }
                  }
                  else {
                      rtype  <- c(rep("E", nr), rep("L", 2*nc))
                  }


                  # ---------------------------------------------
                  # objective function
                  # ---------------------------------------------

                  cobj <- c(rep(0, 2*nc), deltaobj)


                  # ---------------------------------------------
                  # row and column names for the problem object
                  # ---------------------------------------------

                  if (isTRUE(useNames)) {
                      if (is.null(cnames)) {
                          cn <- c(paste("wt",  react_id(model), sep = "_"),
                                  paste("del", react_id(model), sep = "_"),
                                  paste("dM",  react_id(model), sep = "_"),
                                  paste("dP",  react_id(model), sep = "_")
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
                          if (isTRUE(COBRAflag)) {
                              rn <- c(paste("wt",  met_id(model), sep = "_"),
                                      paste("del", met_id(model), sep = "_"),
                                      paste("deltaM", 1:nc, sep = "_"),
                                      paste("deltaP", 1:nc, sep = "_"),
                                      "obj_wt"
                              )
                          }
                          else {
                              rn <- c(paste("del", met_id(model), sep = "_"),
                                      paste("deltaM", 1:nc, sep = "_"),
                                      paste("deltaP", 1:nc, sep = "_")
                              )
                          }
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
                              paste(toupper(subalg), mod_id(model), sep = "_"),
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
                                            sbalg      = subalg,
                                            pType      = "lp",
                                            scaling    = scaling,
                                            fi         = fi,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = LHS,
                                            ub         = upper,
                                            lb         = lower,
                                            obj        = cobj,
                                            rlb        = rlower,
                                            rub        = rupper,
                                            rtype      = rtype,
                                            lpdir      = "min",
                                            ctype      = NULL,
                                            cnames     = colNames,
                                            rnames     = rowNames,
                                            pname      = probName,
                                            algPar     = list("wtflux" = wtflux,
                                                              "wtobj"  = wtobj),
                                            ...)

                  if (!is.null(writeProbToFileName)) {
                      writeProb(problem(.Object),
                                fname = as.character(writeProbToFileName))
                  }


#                  # ---------------------------------------------
#                  # build problem object
#                  # ---------------------------------------------
#
#                  lp <- optObj(solver = solver, method = method)
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
#                             ub    = upper,
#                             lb    = lower,
#                             obj   = cobj,
#                             rlb   = rlower,
#                             rub   = rupper,
#                             rtype = rtype,
#                             lpdir = "min"
#                  )
#                  
#                  if (!is.null(scaling)) {
#                      scaleProb(lp, scaling)
#                  }
#
#                  .Object@problem   <- lp
#                  .Object@algorithm <- subalg
#                  .Object@nr        <- as.integer(nRows)
#                  .Object@nc        <- as.integer(nCols)
#                  .Object@fldind    <- as.integer(fi)
#                  validObject(.Object)
                  
              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
