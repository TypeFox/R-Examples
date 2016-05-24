#  settings.R
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



SYBIL_SETTINGS <- function(parm, value, ...) {

    if ( (missing(parm)) && (missing(value)) ) {
       return(.SYBILenv$settings)
    }

    if (missing(value)) {
        if (!parm %in% names(.SYBILenv$settings)) {
            stop("unknown parameter ", sQuote(parm))
        }
        return(.SYBILenv$settings[[parm]])
    }
    
    if ( (length(parm) != 1) ||
         ( (length(value) != 1) && (! (parm == "SOLVER_CTRL_PARM") ) ) ) {
        stop("arguments 'parm' and 'value' must have a length of 1")
    }
    
    switch(parm,
    
        "SOLVER" = {

            chmet <- checkDefaultMethod(solver = value,
                                        method = "",
                                        probType = "", ...)
            .SYBILenv$settings[["SOLVER"]]           <- chmet$sol
            .SYBILenv$settings[["METHOD"]]           <- chmet$met
            .SYBILenv$settings[["SOLVER_CTRL_PARM"]] <- chmet$parm
        },
    
        "METHOD" = {
            chmet <- checkDefaultMethod(solver = SYBIL_SETTINGS("SOLVER"),
                                        method = value,
                                        probType = "", ...)
            .SYBILenv$settings[["SOLVER"]]           <- chmet$sol
            .SYBILenv$settings[["METHOD"]]           <- chmet$met
            .SYBILenv$settings[["SOLVER_CTRL_PARM"]] <- chmet$parm
        },
    
        "TOLERANCE" = {
            .SYBILenv$settings[["TOLERANCE"]] <- as.numeric(value)
        },
    
        "MAXIMUM" = {
            .SYBILenv$settings[["MAXIMUM"]] <- as.numeric(value)
        },
    
        "ALGORITHM" = {
#            if ( (value == "FBA")            ||
#                 (value == "linearMOMA")     ||
#                 (value == "linearMOMA_COBRA") ) {
#                .SYBILenv$settings[["ALGORITHM"]] <- as.character(value)
#            }
#            else {
#                stop("ALGORITHM can be either 'FBA', ",
#                     "'linearMOMA' or 'linearMOMA_COBRA'")
#            }
            .SYBILenv$settings[["ALGORITHM"]] <- as.character(value)
        },
    
        "OPT_DIRECTION" = {
            if ( (value == "max") || (value == "min") ) {
                .SYBILenv$settings[["OPT_DIRECTION"]] <- as.character(value)
            }
            else {
                stop("OPT_DIRECTION can be either 'max' or 'min'")
            }
        },
    
        "USE_NAMES" = {
            .SYBILenv$settings[["USE_NAMES"]] <- as.logical(value)
        },

        "PATH_TO_MODEL" = {
            if (file.exists(value)) {
                .SYBILenv$settings[["PATH_TO_MODEL"]] <- as.character(value)
            }
            else {
                stop("directory ", sQuote(value), " does not exist")
            }
        },
    
        "SOLVER_CTRL_PARM" = {
            if ( (is.data.frame(value)) || (is.list(value)) ) {
                if ("NA" %in% names(SYBIL_SETTINGS("SOLVER_CTRL_PARM"))) {
                    .SYBILenv$settings[["SOLVER_CTRL_PARM"]] <- value
                }
                else {
                    pn <- names(value)
                    for (i in seq(along = value)) {
                        .SYBILenv$settings[["SOLVER_CTRL_PARM"]][[pn[i]]] <- value[[pn[i]]]
                    }
                }
            }
            else {
                stop("SOLVER_CTRL_PARM must be data.frame or list")
            }
        },

        {
            stop("unknown parameter: ", sQuote(parm))
        }
    )
}

