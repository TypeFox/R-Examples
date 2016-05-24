#  getsybilenv.R
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


getsybilenv <- function(part) {


    # ------------------------------------------------------------------------ #

    printsolvers <- function() {
        cat("\n# --------------------------------------------------------- #\n")
        cat("solver packages:\n")
        cat(paste(.SYBILenv[["solvers"]]), sep = ", ")
        cat("\n\n")
    }


    # ------------------------------------------------------------------------ #

    printmethods <- function() {
        cat("\n# --------------------------------------------------------- #\n")
        cat("methods included in the solver packages:\n")
        slv <- names(.SYBILenv[["solverMethods"]])
        for(i in seq(along = slv)) {
            cat(slv[i], ":\n", sep = "")
            cat(paste(.SYBILenv[["solverMethods"]][[slv[i]]]), sep = ", ")
            cat("\n\n")
        }
        cat("\n")
    }


    # ------------------------------------------------------------------------ #

    printptype <- function() {
        cat("\n# --------------------------------------------------------- #\n")
        cat("methods used for problem types:\n")
        ptype <- names(.SYBILenv[["ptype"]])
        for(i in seq(along = ptype)) {
            cat(ptype[i], ":\n", sep = "")
            slv <- names(.SYBILenv[["ptype"]][[ptype[i]]])
            for(j in seq(along = slv)) {
                cat("    ", slv[j], ":\n    ", sep = "")
                cat(paste(.SYBILenv[["ptype"]][[ptype[i]]][[slv[j]]]), sep = ", ")
                cat("\n\n")
            }
        }
        cat("\n")
    }


    # ------------------------------------------------------------------------ #

    printpurpose <- function() {
        cat("\n# --------------------------------------------------------- #\n")
        cat("algorithms for this purpose:\n")
        fkt <- names(.SYBILenv[["algorithm"]])
        for(i in seq(along = fkt)) {
            cat(fkt[i], ":\n", sep = "")
            cat(paste(.SYBILenv[["algorithm"]][[fkt[i]]]), sep = ", ")
            cat("\n\n")
        }
        cat("\n")
    }


    # ------------------------------------------------------------------------ #

    if (missing(part)) {
        printsolvers()
        printmethods()
        printptype()
        printpurpose()
        #print(.SYBILenv[["solverCtrlParm"]])
    }
    else {
        cmd <- paste("print", part, "()", sep = "")
        err <- tryCatch(eval(parse(text = cmd)), error = function(e) e)
        if (is(err, "simpleError")) {
            stop(sQuote(part), " is not available in the sybil environment")
        }
    }
    
    return(invisible(NULL))

}


