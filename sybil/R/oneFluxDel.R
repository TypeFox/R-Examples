#  oneFluxDel.R
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
# Function: oneFluxDel
#
# This function performs a "gene deletion analysis".
# In each iteration one gene is switched of (vi = 0)
# and the objective function will be computed.
#
# The function oneGeneDel() is inspired by the function
# singleRxnDeletion() contained in the COBRA Toolbox.


oneFluxDel <- function(model, react = c(1:react_num(model)),
                       lb = rep(0, length(react)),
                       ub = rep(0, length(react)),
                       checkOptSolObj = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    creact <- checkReactId(model, react)
    if (!is(creact, "reactId")) {
        stop("check argument react")
    }

    creact <- sort(react_pos(creact))

    sol <- optimizer(model = model, lb = lb, ub = ub,
                     react = as.list(creact), ...)

    # ------------------------------------------------------------------------ #

    optsol <- new("optsol_fluxdel")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt
    
    chlb(optsol) <- as.numeric(lb)
    chub(optsol) <- as.numeric(ub)
    dels(optsol) <- matrix(react_id(model)[creact], ncol = 1)

    if (isTRUE(checkOptSolObj)) {
        checkOptSol(optsol, onlywarn = TRUE)
    }

    return(optsol)

}

