#  printLogComment.R
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
# Function: .printLogComment
#
#
#
#


.printLogComment <- function(reason) {

    # log comment will look like this:

    # -------------------------------------------------------- #
    ## log <reason>: Day Mon  D hh:mm:ss Year
    # -------------------------------------------------------- #

    lc1 <- paste("#", paste(rep("-", 56), collapse = ""), "#", collapse = " ")
    lc2 <- paste("## log", as.character(reason), date())

    lc  <- paste(lc1, lc2, lc1, "\n", sep = "\n")

    return(lc)

}
