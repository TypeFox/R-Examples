#  optObj_basicfunc.R
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

# print a warning
wrong_type_msg <- function(lp) {
    warning(paste("Slot oobj of", lp@solver,
                  "is not a valid pointer to a valid solver."
            )
    )
}


wrong_solver_msg <- function(lp, method, printOut = TRUE) {
    if (isTRUE(printOut)) {
        warning(paste("Requested method", method,
                      "is not available for solver", lp@solver
                )
        )
    }
}


#------------------------------------------------------------------------------#

# which optimizations did not end successful
checkSolStat <- function(stat, solver = SYBIL_SETTINGS("SOLVER")) {

    out <- FALSE
    switch(solver,
        "glpkAPI" = {
            out <- which(stat != 5)
        },
        "clpAPI" = {
            out <- which(stat != 0)
        },
        "lpSolveAPI" = {
            out <- which(stat != 0)
        },
        "cplexAPI" = {
            # CPLEX: 101 optimal integer solution
            # CPLEX: 102 Optimal solution with the tolerance defined by epgap or epagap
            out <- which(stat != 1 & stat != 101 & stat != 102)
        },
        {
            cmd <- paste(solver,"::checkSolutionStatus(stat)", sep = "")
            out <- tryCatch(eval(parse(text = cmd)), error = function(e) NA)    
            #warning("not a valid solver")
        }
    )
    return(out)
}


#------------------------------------------------------------------------------#

getMeanReturn <- function(code, solver = SYBIL_SETTINGS("SOLVER")) {

    out <- FALSE
    switch(solver,
        "glpkAPI" = {
            out <- glpkAPI::return_codeGLPK(code)
        },
        "clpAPI" = {
            out <- clpAPI::return_codeCLP(code)
        },
        "lpSolveAPI" = {
            out <- return_codeLPSOLVE(code)
        },
        "cplexAPI" = {
            out <- cplexAPI::return_codeCPLEX(code)
        },
        {
            cmd <- paste(solver,"::getReturnString(code)", sep = "")
            out <- tryCatch(eval(parse(text = cmd)),
                            error = function(e) as.character(NA))    
            #warning("not a valid solver")
        }
    )
    return(out)
}


#------------------------------------------------------------------------------#

getMeanStatus <- function(code,
                          solver = SYBIL_SETTINGS("SOLVER"), env = NULL) {
    out <- FALSE
    switch(solver,
        "glpkAPI" = {
            out <- glpkAPI::status_codeGLPK(code)
        },
        "clpAPI" = {
            out <- clpAPI::status_codeCLP(code)
        },
        "lpSolveAPI" = {
            #out <- "see return code"
            out <- getMeanReturn(code, solver)
        },
        "cplexAPI" = {
            if (is.null(env)) {
                cenv <- cplexAPI::openEnvCPLEX()
            }
            else {
                cenv <- env
            }
            out <- cplexAPI::status_codeCPLEX(cenv, code)
            if (is.null(env)) {
                cplexAPI::closeEnvCPLEX(cenv)
            }
            rm(cenv)
        },
        {
            cmd <- paste(solver,"::getStatusString(code)", sep = "")
            out <- tryCatch(eval(parse(text = cmd)),
                            error = function(e) as.character(NA))    
            #warning("not a valid solver")
        }
    )
    return(out)
}


#TIMEOUT <- "timeout"
#INFINITE <- "infinite"
