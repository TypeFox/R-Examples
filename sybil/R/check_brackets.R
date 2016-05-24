#  check_brackets.R
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
# Function: .check_brackets
#
#
#


# check, whether we have the same number of '(' and ')'.
.check_brackets <- function(rule) {

    i  <- 1
    nb <- 0
    nc <- length(rule)

    while ( (nb >= 0) && (i <= nc) ) {
        switch(rule[i],
            "(" = { nb <- nb + 1 },
            ")" = { nb <- nb - 1 }
        )
        i <- i + 1
    }

    return( ifelse(nb == 0, TRUE, FALSE) )

}
