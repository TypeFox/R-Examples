#  doInRound.R
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
# Function: .doInRound
#
#
#


.doInRound <- function(indRounds, maxRounds) {

    if (maxRounds < 1) {
        stop("Argument 'maxRounds' must be > 0")
    }

    if (is.null(indRounds)) {
        runPP <- c(1:maxRounds)
    }
    else {
        if (is(indRounds, "character")) {

            # make indRounds integer values
            DIR <- as.integer(indRounds)

            # first value is treated as offset, if indRounds
            # contains two elements (more than two elements are ignored)
            if (length(DIR) > 1) {
                offs <- DIR[1]
                DIR  <- ifelse(DIR[2] < 1, 1L, abs(DIR[2]))
            }
            else {
                offs <- 0L
                DIR  <- ifelse(DIR[1] < 1, 1L, abs(DIR[1]))
            }

            # when we will run pre/post processing
            runPP <- seq(from = (DIR+offs), to = maxRounds, by = DIR)
            runPP <- runPP[runPP > 0]

        }
        else if (is(indRounds, "numeric")) {
            runPP <- sort(as.integer(indRounds[indRounds > 0]))
        }
        else {
            warning("Argument 'indRounds' must be numeric or character")
            runPP <- c(1:maxRounds)
        }
    }

    return(runPP)

}

