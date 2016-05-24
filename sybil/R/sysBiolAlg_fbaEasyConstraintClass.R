#  sysBiolAlg_fbaEasyConstraintClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Copyright (C) 2014-2015 Claus Jonathan Fritzemeier, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: clausjonathan.fritzemeier@hhu.de
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
#          definition of the class sysBiolAlg_fbaEasyConstraint                #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_fbaEasyConstraint",
			representation(
				easyConstraint = "list"
			),
			contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_fbaEasyConstraint
setMethod(f = "initialize",
          signature = "sysBiolAlg_fbaEasyConstraint",
          definition = function(.Object,
                                model,
                                lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                easyConstraint = NULL,
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
                  
                  mat <- S(model)
                  rtype <- rep("E", nRows)
                  rlb <- rep(0, nRows)
                  rub <- rep(0, nRows)
                  
                  #add easyConstraints:
                  if(!is.null(easyConstraint)){
                  	if(		length(easyConstraint$react) != length(easyConstraint$x)
                  		| 	length(easyConstraint$react) != length(easyConstraint$rtype)
                  		){
                  		stop("easyConstraint elements have to have equal lengths")
                  	}
                  	stopifnot(is.list(easyConstraint$react))
                  	stopifnot(is.list(easyConstraint$x))
                  	stopifnot(all(easyConstraint$rtype %in% c("F", "L", "U", "D", "E")))
                  	
                  	# setting and checking rlb
                  	if(is.null(easyConstraint$lb)){
                  		rlb <- c(rlb, rep(0, length(easyConstraint$react)))
                  	}else{
                  		if(length(easyConstraint$react) != length(easyConstraint$lb)){
                  			stop("easyConstraint$lb length has to match length of react argument")
                  		}else{
                  			stopifnot(is.numeric(easyConstraint$lb))
                  			rlb <- c(rlb, easyConstraint$lb)
                  		}
                  	}
                  	
                  	# setting and checking rub
                  	if(is.null(easyConstraint$ub)){
                  		rub <- c(rub, rep(0, length(easyConstraint$react)))
                  	}else{
                  		if(length(easyConstraint$react) != length(easyConstraint$ub)){
                  			stop("easyConstraint$ub length has to match length of react argument")
                  		}else{
                  			stopifnot(is.numeric(easyConstraint$ub))
                  			rub <- c(rub, easyConstraint$ub)
                  		}
                  	}
                  	
                  	m <- Matrix(0, ncol=nCols, nrow=length(easyConstraint$react))
                  	
                  	for(i in 1:length(easyConstraint$react)){
                  		m[i, easyConstraint$react[[i]]] <- easyConstraint$x[[i]]
                  	}
                  	
                  	
                  	mat <- rbind2(mat, m)
                  	rtype <- c(rtype, easyConstraint$rtype)
                  	nRows <- nRows + length(easyConstraint$react)
                  	if(!is.null(rowNames)){
                  		rowNames <- c(rowNames, paste0("easyConstraint", 1:length(easyConstraint$react)))
                  	}
                  	
                  }
                  
                  # generate problem object
                  .Object <- callNextMethod(.Object,
                                            sbalg      = "fba",
                                            pType      = "lp",
                                            scaling    = scaling,
                                            fi         = 1:nCols,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = mat,
                                            ub         = uppbnd(model),
                                            lb         = lowbnd(model),
                                            obj        = obj_coef(model),
                                            rlb        = rlb,
                                            rtype      = rtype,
                                            lpdir      = lpdir,
                                            rub        = rub,
                                            ctype      = NULL,
                                            cnames     = colNames,
                                            rnames     = rowNames,
                                            pname      = probName,
                                            ...)
					.Object@easyConstraint <- easyConstraint
                  if (!is.null(writeProbToFileName)) {
                      writeProb(problem(.Object),
                                fname = as.character(writeProbToFileName))
                  }
              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
