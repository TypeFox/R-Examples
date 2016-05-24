#  createReactionString.R
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
# Function: .createReactionString
#
#
#

.createReactionString <- function(model,
                                  makeClosedNetwork,
                                  entrydelim = ", ",
                                  extMetFlag = "b") {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    # required data structures
    equat  <- vector(mode = "character", length = react_num(model))
    compat <- vector(mode = "character", length = react_num(model))

    revers <- ifelse(react_rev(model), "Reversible", "Irreversible")
    arrow  <- ifelse(react_rev(model), "<==>", "-->")


    # remove compartment flag if existing
    metab <- sub("\\[\\w+\\]$", "", met_id(model))
    metcp <- sub(".+(\\[\\w+\\])$", "\\1", met_id(model))

    for (j in 1:react_num(model)) {
        column <- S(model)[,j]

        # row indices
        constr_ind <- which(column != 0)
        # stoichiometric coefficients
        stcoef     <- column[constr_ind]

        # check if reaction is empty
        if (length(constr_ind) > 0) {

            comp <- unique(mod_compart(model)[met_comp(model)[constr_ind]])

            # reaction involves more than one compartment
            if (length(comp) > 1) {
                if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
                    metr      <- paste(met_id(model)[constr_ind],
                                       "[",
                                       mod_compart(model)[met_comp(model)[constr_ind]],
                                       "]",
                                       sep = "")
                }
                else {
                    metr      <- met_id(model)[constr_ind]
                }
                compat[j] <- paste(comp, collapse = entrydelim)
                compflag  <- ""
            } else {
                # Check if the current reaction is an exchange reaction.
                # In order to build a closed network, we need to add a 'boundary'
                # metabolite [b].
                if ( (isTRUE(makeClosedNetwork)) && (length(constr_ind) == 1) ) {
                    if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
                        metIN     <- paste(met_id(model)[constr_ind],
                                           "[",
                                           mod_compart(model)[met_comp(model)[constr_ind[1]]],
                                           "]",
                                           sep = "")
                    }
                    else {
                        metIN     <- met_id(model)[constr_ind]
                    }
                    metr       <- c(metIN,
                                    paste(metab[constr_ind],
                                          "[", extMetFlag, "]",sep = ""))
                    constr_ind <- c(constr_ind, constr_ind)
                    stcoef     <- c(stcoef, (stcoef * -1))
                    compat[j]  <- comp
                    compflag   <- ""
                }
                else {
                    metr      <- metab[constr_ind]
                    compat[j] <- comp
                    # if yes, the metabolite id does not contain the metabolite
                    # compartment
                    if (metcp[constr_ind[1]] == metab[constr_ind[1]]) {
                        compflag  <- paste("[",
                                           mod_compart(model)[met_comp(model)[constr_ind[1]]],
                                           "]",
                                           sep = "")
                    }
                    else {
                        compflag  <- metcp[constr_ind[1]]
                    }
                }
            }

            educt <- vector(mode = "list")
            product <- vector(mode = "list")
    
            for (i in seq(along = constr_ind)) {
                if (stcoef[i] > 0) {
                    stoich <- ifelse(stcoef[i] != 1,
                                     paste("(", stcoef[i], ") ", sep = ""),
                                     "")
                    product[metr[i]] <- paste(stoich, metr[i], sep = "")
                }
                else {
                    stoich <- ifelse(stcoef[i] != -1,
                                     paste("(", (stcoef[i] * -1), ") ", sep = ""),
                                     "")
                    educt[metr[i]] <- paste(stoich, metr[i], sep = "")
                }
            }
    
    
            equattmp <- paste(paste(educt, collapse = " + "),
                              arrow[j],
                              paste(product, collapse = " + "))
    
            if (compflag == "") {
                equat[j] <- equattmp
            }
            else {
                equat[j] <- paste(compflag, equattmp, sep = " : ")
            }
    
        }

    }

    return(list(equat  = equat,
                compat = compat,
                revers = revers,
                metab  = metab))

}
