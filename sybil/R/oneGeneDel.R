#  oneGeneDel.R
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
# Function: oneGeneDel
#
# This function performs a "gene deletion analysis".
# In each iteration all reactions corresponding to a
# gene were switched of (vi = 0) and the objective
# function will be computed.
#
# The function oneGeneDel() is inspired by the function
# singleGeneDeletion() contained in the COBRA Toolbox.


oneGeneDel <- function(model, geneList,
                       lb = rep(0, length(geneList)),
                       ub = rep(0, length(geneList)),
                       checkOptSolObj = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    if (missing(geneList)) {
        if (length(allGenes(model)) < 1) {
            stop("Argument 'geneList' must contain at least one gene!")
        }
        else {
            geneList <- c(1:length(allGenes(model)))
        }
    }
 
    # translate the gene List in indices of allGenes(model)
    if (is(geneList, "character")) {
       geneList <- match(geneList, allGenes(model))
       if (any(is.na(geneList))) {
           stop("check genelist!")
       }
    }

    #geneList <- sort(geneList)


    # ------------------------------------------------------------------------ #

    fd <- .generateFluxdels(model, geneList)

    sol <- optimizer(model = model,
                     react = fd[["react"]], lb = lb, ub = ub, ...)


    # ------------------------------------------------------------------------ #

    optsol <- new("optsol_genedel")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt
    
    chlb(optsol)      <- as.numeric(lb)
    chub(optsol)      <- as.numeric(ub)
    dels(optsol)      <- matrix(allGenes(model)[geneList], ncol = 1)
    fluxdels(optsol)  <- fd[["fd"]]
    hasEffect(optsol) <- fd[["heff"]]

    if (isTRUE(checkOptSolObj)) {
        checkOptSol(optsol, onlywarn = TRUE)
    }

    return(optsol)

}



