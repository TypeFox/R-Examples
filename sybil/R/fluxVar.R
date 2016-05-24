#  fluxVar.R
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
# Function: fluxVar
#
# Performs a flux vaiability Analysis
# 
# The function fluxVar() is inspired by the function
# fluxVariability() contained in the COBRA Toolbox.
# The algorithm is the same.


fluxVar <- function(model, react = c(1:react_num(model)), exex = FALSE, ...) {

    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }

    # remove exchange reactions from analysis
    if (isTRUE(exex)) {
        exchReact <- findExchReact(model)
        ex <- react_pos(exchReact)
        intReact <- 1:react_num(model)
        intReact <- intReact[-ex]
        
        if (length(intReact) < 1) {
            stop("model contains no internal reactions!")
        }
    }
    else {
        intReact <- react
    }

    creact <- checkReactId(model, intReact)
    
    if (!is(creact, "reactId")) {
        stop("check argument react")
    }


#------------------------------------------------------------------------------#
#                               optimizations                                  #
#------------------------------------------------------------------------------#

    sol <- optimizer(model,
                     react = as.list(c(react_pos(creact), react_pos(creact))),
                     obj_coef = rep(1, (2*length(creact))),
                     lpdir = c(rep("min", length(creact)),
                               rep("max", length(creact))),
                     algorithm = "fv", ...)
    

#------------------------------------------------------------------------------#
#                             save the results                                 #
#------------------------------------------------------------------------------#

    optsol <- new("optsol_fluxVar")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt

    react(optsol) <- creact

    return(optsol)

}


 
