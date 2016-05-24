#  phpp.R
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
# Function: phpp
#
#
# The function phpp() is inspired by the function
# phenotypePhasePlane() contained in the COBRA Toolbox.
# The algorithm is the same.


phpp <- function(model, ctrlreact, rng = c(0, 0, 20, 20),
                 numP = 50, setToZero = TRUE, redCosts = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    
    if (length(ctrlreact) != 2) {
        stop("Please enter two control reactions.")
    }
    
    tmp <- checkReactId(model, ctrlreact)
    
    if (!is(tmp, "reactId")) {
        stop("Check control reaction!")
    }
    
    ctrlr <- react_pos(tmp)

    if (length(rng) != 4) {
        stop("rng must have length 4")
    }

    rngA <- seq(rng[1], rng[3], length.out = numP)
    rngB <- seq(rng[2], rng[4], length.out = numP)
    ctrlfl <- as.matrix(expand.grid(rngA, rngB))
    colnames(ctrlfl) <- react_id(tmp)

    num_of_prob <- numP ** 2

    # reduced costs
    ca <- match.call()
    if ("poCmd" %in% names(ca)) {
        rccmd <- ca[["poCmd"]]
    }
    else {
        if (isTRUE(redCosts)) {
            rccmd <- list("getRedCosts")
        }
        else {
            rccmd <- list(NA)
        }
    }

#------------------------------------------------------------------------------#
#                                optimizations                                 #
#------------------------------------------------------------------------------#

    sol <- optimizer(model,
                     react = rep(list(ctrlr), num_of_prob),
                     lb = -1*abs(ctrlfl),
                     ub = -1*abs(ctrlfl),
                     algorithm = "fba",
                     setToZero = setToZero,
                     poCmd = rccmd, ...)


#------------------------------------------------------------------------------#
#                             save the results                                 #
#------------------------------------------------------------------------------#

    optsol <- new("optsol_phpp")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt

    ctrlr(optsol)  <- tmp
    ctrlfl(optsol) <- as.matrix(ctrlfl)

    if (isTRUE(redCosts)) {
        rc  <- pa(postProc(optsol))
        src <- matrix(unlist(rc),
                      nrow = length(rc), byrow = TRUE)[, ctrlr, drop = FALSE]
        if (isTRUE(setToZero)) {
            tz <- checkStat(optsol)
            src[tz, ] <- 0
        }
        colnames(src) <- react_id(tmp)
        optsol@redCosts <- src
        np <- ppProc(list(NA))
        pa(np) <- list(NA)
        postProc(optsol) <- np
    }
    else {
        optsol@redCosts <- as.matrix(NA)
    }

    return(optsol)

}



