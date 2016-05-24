#  changeBounds.R
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
# Function: changeBounds
#
#
#
#

changeBounds <- function(model, react, lb = NULL, ub = NULL) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
  
    checkedIds <- checkReactId(model, react)
    if (!is(checkedIds, "reactId")) {
        stop("argument react is wrong")
    }

    if ( (is.null(lb)) && (is.null(ub)) ) {
        lowbnd(model)[react_pos(checkedIds)] <- rep(0, length(checkedIds))
        uppbnd(model)[react_pos(checkedIds)] <- rep(0, length(checkedIds))
    }
    else {
        # set upper bound
        if (!is.null(ub)) {
            stopifnot(is(ub, "numeric"))
            if (length(ub) == 1) {
                ubnd <- rep(ub, length(checkedIds))
            }
            else {
                stopifnot(length(ub) == length(checkedIds))
                ubnd <- ub
            }
            uppbnd(model)[react_pos(checkedIds)] <- ubnd
        }

        if (!is.null(lb)) {
            stopifnot(is(lb, "numeric"))
            if (length(lb) == 1) {
                lbnd <- rep(lb, length(checkedIds))
            }
            else {
                stopifnot(length(lb) == length(checkedIds))
                lbnd <- lb
            }
            lowbnd(model)[react_pos(checkedIds)] <- lbnd
        }
    }

    return(model)

}

