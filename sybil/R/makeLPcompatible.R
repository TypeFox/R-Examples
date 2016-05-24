#  makeLPcompatible.R
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
# Function: .makeLPcompatible
#
# generates names, compatible to CPLEX LP file format:
# a name must not start with a period or digit and must not contain
# square brackets

.makeLPcompatible <- function(name, prefix, sep = "") {

    nameLP <- sub("[", "(", name,   fixed = TRUE)
    nameLP <- sub("]", ")", nameLP, fixed = TRUE)
    nameLP <- sub("^([0-9\\.])", paste(prefix, "\\1", sep = sep),
                  nameLP, perl = TRUE)

    #nameLP <- paste(prefix, nameLP, sep = sep) # prefix = "M_" or "R_"

    return(nameLP)

}
