#  validreactId.R
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
# Function: .validreactId
#
# Validity checking of an object of class reactId
#


.validreactId <- function(object) {

  if (!is(object, "reactId")) {
      return("needs an object of class reactId!")
  }

  if (length(mod_id(object)) != 1) {
      return("slot mod_id must be of length 1")
  }

  if (length(react_pos(object)) != length(object)) {
      return(paste("slot react_pos must be of length", length(object)))
  }

  if (length(react_id(object)) != length(object)) {
      return(paste("slot react_id must be of length", length(object)))
  }

  return(TRUE)
}

