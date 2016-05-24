#  editEnvir.R
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
# Function: editEnvir
#
#
#

editEnvir <- function(model, newKey = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("model must be of class modelorg")
    }

    ex <- findExchReact(model)
    exid <- react_pos(ex)

    exfr <- data.frame("reaction_id"   = react_id(model)[exid],
                       "lower_bound"   = lowbnd(model)[exid],
                       "upper_bound"   = uppbnd(model)[exid],
                       "reaction_name" = react_name(model)[exid])

    exfr <- edit(exfr, ...)

    lowbnd(model)[exid] <- exfr[["lower_bound"]]
    uppbnd(model)[exid] <- exfr[["upper_bound"]]
    
    if (isTRUE(newKey)) {
        mod_key(model) <- .generateModKey()
    }

    validObject(model)

    return(model)

}
