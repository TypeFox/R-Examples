#  changeObjFunc.R
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
# Function: changeObjFunc
#
#
#
#

changeObjFunc <- function(model, react, obj_coef = rep(1, length(react))) {

  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }

  if (length(react) != length(obj_coef)) {
      stop("react and obj_coef must have the same length!")
  }

  checkedIds <- checkReactId(model, react)
  if (!is(checkedIds, "reactId")) {
      stop("Check your reaction Id's")
  }
 
  # set all objective coefficients to zero
  obj_coef(model) <- numeric(react_num(model))

  obj_coef(model)[react_pos(checkedIds)] <- obj_coef

  return(model)

}

