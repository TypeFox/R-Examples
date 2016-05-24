#  optObj_lpSolveAPIcompat.R
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
# This is for compatibility with lpSolveAPI.
# lpSolveAPI only has a return code, which also acts as status code.

return_codeLPSOLVE <- function(code) {
    if (code == 0)       { return( "optimal solution found" ) }
    else if (code == 1)  { return( "the model is sub-optimal" ) }
    else if (code == 2)  { return( "the model is infeasible" ) }
    else if (code == 3)  { return( "the model is unbounded" ) }
    else if (code == 4)  { return( "the model is degenerate" ) }
    else if (code == 5)  { return( "numerical failure encountered" ) }
    else if (code == 6)  { return( "process aborted" ) }
    else if (code == 7)  { return( "timeout" ) }
    else if (code == 9)  { return( "the model was solved by presolve" ) }
    else if (code == 10) { return( "the branch and bound routine failed" ) }
    else if (code == 11) { return( paste("the branch and bound was stopped",
                                         "because of a break-at-first",
                                         "or break-at-value"
                                   )
                           )
    }
    else if (code == 12) { return( paste("a feasible branch and bound",
                                         "solution was found"
                                   )
                           )
    }
    else if (code == 13) { return( paste("no feasible branch and bound",
                                         "solution was found"
                                   )
                           )
    }
    else { return(paste("Failed to obtain solution",
                        "unknown error code:", code
                  )
           )
    }
}


#------------------------------------------------------------------------------#

loadMatrixPerColumnLPSOLVE <- function(lpmod, constMat) {

    stopifnot(is(constMat, "Matrix"))
    
    x <- constMat@x
    p <- constMat@p + 1
    i <- constMat@i + 1

    k <- 1
    while (k <= ncol(constMat)) {
        lpSolveAPI::set.column(lpmod,
                               column  = k,
                               x       = x[(p[k]):(p[k+1]-1)],
                               indices = i[(p[k]):(p[k+1]-1)])
        k <- k + 1
    }

}


#------------------------------------------------------------------------------#

# finalizeLpSolveProb <- function(lprec) {
# 
#     if (is(lprec, "lpExtPtr")) {
#         #tryCatch(#.Call("RlpSolve_delete_lp", lprec, PACKAGE = "lpSolveAPI"),
#         #    lpSolveAPI::delete.lp(lprec),
#         #    warning = function(e) FALSE)
#         capture.output(lpSolveAPI::delete.lp(lprec), file = "/dev/null")
#     }
#     invisible(lprec)
# }


