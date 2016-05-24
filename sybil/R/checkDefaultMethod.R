#  checkDefaultMethod.R
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


checkDefaultMethod <- function(solver, method, probType, loadPackage = TRUE) {

    stopifnot(is(solver, "character"),
              is(method, "character"),
              is(probType, "character"))

#    # -------------------------------------------------------------- #
#    # load solver package
#
#    if (isTRUE(loadPackage)) {
#        checkPackage <- require(solver, character.only = TRUE)
#    
#        if(!isTRUE(checkPackage)) {
#            stop("package ", sQuote(solver), " not found")
#        }
#    }

    # -------------------------------------------------------------- #
    # validate solver

    val_solver_ind <- match(solver, .SYBILenv$solvers)

    if (is.na(val_solver_ind)) {
        cmd <- paste(solver, ".onAttach()", sep = ":::")
        slv <- tryCatch(eval(parse(text = cmd)), error = function(e) e)
        if (isTRUE(slv)) {
            val_solver_ind <- match(solver, .SYBILenv$solvers)
        }

        if (is.na(val_solver_ind)) {
            val_solver_ind <- 1L
            warning("solver ", sQuote(solver),
                    " not found, using default: ",
                    sQuote(.SYBILenv$solvers[val_solver_ind]))
        }
    }

    val_solver <- .SYBILenv$solvers[val_solver_ind]


    # -------------------------------------------------------------- #
    # validate method

    val_method_ind <- match(method, .SYBILenv$solverMethods[[val_solver]])

    if (is.na(val_method_ind)) {
        val_method <- .SYBILenv$solverMethods[[val_solver]][1]
    }
    else {
        val_method <- .SYBILenv$solverMethods[[val_solver]][val_method_ind]
    }


    # -------------------------------------------------------------- #
    # validate method with problem type

    if (probType %in% names(.SYBILenv$ptype)) {
        if (val_solver %in% names(.SYBILenv$ptype[[probType]])) {
            meth_tmp <- match(val_method, .SYBILenv$ptype[[probType]][[val_solver]])
            if (is.na(meth_tmp)) {
                val_method <- .SYBILenv$ptype[[probType]][[val_solver]][1]
            }
        }
        else {
            stop("solver ", sQuote(val_solver),
                 " can not handle problems of type ", sQuote(probType))
        }
    }


    # -------------------------------------------------------------- #
    # solver parameters

    ctrl_parm <- .SYBILenv$solverCtrlParm[[val_solver]][[val_method]]

    if (is.null(ctrl_parm)) {
        ctrl_parm <- as.data.frame(NA)
    }


    # -------------------------------------------------------------- #
    # load solver package

    if (isTRUE(loadPackage)) {
        checkPackage <- require(val_solver, character.only = TRUE)
    
        if(!isTRUE(checkPackage)) {
            stop("package ", sQuote(val_solver), " not found")
        }
    }

    # -------------------------------------------------------------- #
    # done

    return(list(sol  = as.character(val_solver),
                met  = as.character(val_method),
                parm = ctrl_parm
           )
    )

}

