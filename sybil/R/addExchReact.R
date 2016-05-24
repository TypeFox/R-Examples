#  addExchReact.R
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
# Function: addExchReact
#
#
# The function addExchReact() is inspired by the function
# addExchangeRxn() contained in the COBRA Toolbox.
# The algorithm is (more or less) the same.


addExchReact <- function(model, met, lb, ub) {

  
    # ------------------------------------------------------------------------ #
    # check arguments
    # ------------------------------------------------------------------------ #

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg")
    }

    if ( (length(met) < 1) || (met == "") ) {
        stop("at least one metabolite is required")
    }

    if (missing(lb)) {
        Clb <- rep(0, length(met))
    }
    else {
        Clb <- lb
    }

    if (missing(ub)) {
        Cub <- rep(SYBIL_SETTINGS("MAXIMUM"), length(met))
    }
    else {
        Cub <- ub
    }

    if ( (length(met) != length(Clb)) || (length(met) != length(Cub)) ) {
        stop("arguments 'met', 'lb' and 'ub' must have the same length")
    }

    Crev <- rep(FALSE, length(met))
    Crev[Clb < 0] <- TRUE
    
    exRid <- paste("Ex_", met, sep = "")
    
    mod_out <- model
    
    for (i in seq(along = met)) {
        mod_out <- addReact(model = mod_out,
                            id = exRid[i],
                            met = met[i],
                            Scoef = -1,
                            reversible = Crev[i],
                            lb = Clb[i],
                            ub = Cub[i])
    }

    return(mod_out)
}

