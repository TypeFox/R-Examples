#  robAna.R
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
# Function: robAna
#
#
# The function robAna() is inspired by the function
# robustnessAnalysis() contained in the COBRA Toolbox.
# The algorithm is the same.


robAna <- function(model, ctrlreact, rng = NULL,
                   numP = 20, verboseMode = 1, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    
    if (length(ctrlreact) != 1) {
        stop("Please enter exactly one control reaction.")
    }
    
    tmp <- checkReactId(model, ctrlreact)
    
    if (!is(tmp, "reactId")) {
        stop("Check control reaction!")
    }
    
    ctrlr <- react_pos(tmp)


#------------------------------------------------------------------------------#
#                       minimum and maximum solution                           #
#------------------------------------------------------------------------------#

    if (is.null(rng)) {
        suppressMessages(
            fvmm <- fluxVar(model,
                            react = ctrlr,
                            fixObjVal = FALSE,
                            verboseMode = 0, ...)
        )
        if (length(checkStat(fvmm)) > 0) {
            stop("Optimization for min/max solution ended not successfull!")
        }
        else {
            mm <- c(minSol(fvmm, "lp_obj"), maxSol(fvmm, "lp_obj"))
        }
    }
    else {
        stopifnot(length(rng) == 2, is(rng, "numeric"))
        mm <- sort(rng)
    }

    # sequence of numP numbers between lpmin and lpmax,
    # all with the same distance
    ctrlfl <- seq(mm[1], mm[2], length.out = numP)


#------------------------------------------------------------------------------#
#                                optimization                                  #
#------------------------------------------------------------------------------#

    sol <- optimizer(model,
                     react = as.list(rep(ctrlr, numP)),
                     lb = ctrlfl,
                     ub = ctrlfl,
                     algorithm = "fba",
                     verboseMode = verboseMode, ...)


#------------------------------------------------------------------------------#
#                             save the results                                 #
#------------------------------------------------------------------------------#

    optsol <- new("optsol_robAna")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt

    ctrlr(optsol)  <- tmp
    ctrlfl(optsol) <- as.numeric(ctrlfl)

    return(optsol)

}



