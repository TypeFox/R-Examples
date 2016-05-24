#  validmodelorg.R
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
# Function: .validmodelorg
#
# Validity checking of an object of class modelorg.
#
# Returns TRUE if the model is valid, otherwise
# a character String containing a description of
# the error.


.validmodelorg <- function(object) {

    if (!is(object, "modelorg")) {
        return("needs an object of class modelorg!")
    }
    
    if ((length(mod_id(object)) != 1) || (length(mod_name(object)) != 1)) {
        return("mod_id and mod_name must have a length of 1!")
    }
    
    if (length(mod_desc(object)) == 0) {
    
        if (length(mod_name(object)) != 1) {
            return("mod_name must have a length of 1!")
        }
        if (length(mod_id(object)) != 1) {
            return("mod_id must have a length of 1!")
        }
    
    }
    else {
    
        # model describing stuff
        if (length(mod_desc(object)) != 1) {
            return("description must have a length of 1!")
        }
        if (length(mod_name(object)) != 1) {
            return("mod_name must have a length of 1!")
        }
        if (length(mod_id(object)) != 1) {
            return("mod_id must have a length of 1!")
        }
        if (length(mod_compart(object)) != length(unique(mod_compart(object)))) {
            dup <- duplicated(met_id(object))
            return(paste("mod_compart must be unique! Check ", paste(mod_compart(object)[dup], collapse = ", "), ".", sep = ""))
        }
    
        # metabolite stuff
        if (length(met_num(object)) != 1) {
            return("met_num must have a length of 1!")
        }
        met <- met_num(object)
        if (met != length(met_id(object))) {
            return("Wrong number of metabolite id's!")
        }
        if (met != length(unique(met_id(object)))) {
            dup <- duplicated(met_id(object))
            return(paste("met_id must be unique! Check ", paste(met_id(object)[dup], collapse = ", "), ".", sep = ""))
        }
        if (met != length(met_name(object))) {
            return("Wrong number of metabolite names!")
        }
        if (met != length(met_comp(object))) {
            return("Wrong number of metabolite compartments!")
        }
        if (met != length(met_single(object))) {
            return("Wrong length of met_single!")
        }
        if (met != length(met_de(object))) {
            return("Wrong length of met_de!")
        }
    
        # reactions stuff
        if (length(react_num(object)) != 1) {
            return("react_num must have a length of 1!")
        }
        react <- react_num(object)
        if (react != length(react_id(object))) {
            return("Wrong number of reaction id's!")
        }
        if (react != length(unique(react_id(object)))) {
            dup <- duplicated(react_id(object))
            return(paste("reaction_id must be unique! Check ", paste(react_id(object)[dup], collapse = ", "), ".", sep = ""))
        }
        if (react != length(react_name(object))) {
            return("Wrong number of reaction names!")
        }
        if (react != length(react_rev(object))) {
            return("Wrong number of reversibilities!")
        }
        if (react != length(obj_coef(object))) {
            return("Wrong length of lower bounds!")
        }
        if (react != length(lowbnd(object))) {
            return("Wrong length of lower bounds!")
        }
        if (react != length(uppbnd(object))) {
            return("Wrong length of upper bounds!")
        }
        if (react != length(react_single(object))) {
            return("Wrong length of react_single!")
        }
        if (react != length(react_de(object))) {
            return("Wrong length of react_de!")
        }

        # stoichiometric matrix
        if (identical(dim(S(object)), c(met, react)) == FALSE) {
            return("Wrong dimension of S!")
        }
    
        # GPR stuff
        if (length(gprRules(object)) != 0) {

            if (react != length(gprRules(object))) {
                return("Wrong length of gprRules!")
            }
            if (react != length(genes(object))) {
                return("Wrong length of genes!")
            }
            if (react != length(gpr(object))) {
                return("Wrong length of GPR associations!")
            }
            if (react != nrow(subSys(object))) {
                return("Wrong number of sub systems!")
            }

            # It is possible, that the number of unique genes is smaller than
            # the number of reactions. So, I commented the following out.
            ##print(react)
            ##print(length(allGenes(object)))
            #if (react > length(allGenes(object))) {
            #    return("Wrong number of unique genes!")
            #}

            #print(dim(rxnGeneMat(object)))
            #print(c(react, length(allGenes(object))))
            if (identical(dim(rxnGeneMat(object)), c(react, length(allGenes(object)))) == FALSE) {
                return("Wrong dimension of rxnGeneMat!")
            }
        }
    }
    return(TRUE)
}
