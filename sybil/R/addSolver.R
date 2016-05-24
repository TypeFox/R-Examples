#  addSolver.R
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


################################################
# Function: addSolver
#
#
#

addSolver <- function(solver, method, probType) {
  
    stopifnot(inherits(solver, "character"),
              inherits(method, "character"),
              inherits(probType, "list"),
              (length(solver) == 1L),
              length(probType) == length(method))


    # -------------------------------------------------------------------- #
    # solver
    if (solver %in% .SYBILenv$solvers) {
        return(invisible(FALSE))
    }
    else {
        .SYBILenv$solvers <- append(.SYBILenv$solvers, solver)
    }


    # -------------------------------------------------------------------- #
    # method
    .SYBILenv$solverMethods[[solver]] <- method
    .SYBILenv$solverCtrlParm[[solver]] <- vector(length = length(method),
                                                 mode = "list")

    names(.SYBILenv$solverCtrlParm[[solver]]) <- method

    for (i in seq(along = method)) {
        .SYBILenv$solverCtrlParm[[solver]][[i]] <- as.data.frame(NA)
    }


    # -------------------------------------------------------------------- #
    # problem type
    for (i in seq(along = method)) {
        pt <- unique(probType[[i]])
        stopifnot(inherits(pt, "character"))
        prtp <- pt %in% names(.SYBILenv$ptype)
        for (j in seq(along = pt[prtp])) {
            .SYBILenv$ptype[[pt[prtp][j]]][[solver]] <- append(.SYBILenv$ptype[[pt[prtp][j]]][[solver]], method[i])
        }
    }

    return(invisible(TRUE))

}

