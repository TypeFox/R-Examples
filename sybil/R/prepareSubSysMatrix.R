#  prepareSubSysMatrix.R
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
# Function: .prepareSubSysMatrix
#
#
#

.prepareSubSysMatrix <- function(reactSUBS, nreact, entrydelim = NA) {

    if (!is.null(reactSUBS)) {

        reactSsys <- as.character(reactSUBS)

        ssys_tmp <- strsplit(reactSsys, entrydelim, fixed = TRUE)

        sys_names <- unique(unlist(ssys_tmp))
        ssys <- Matrix::Matrix(FALSE,
                               nrow = length(reactSUBS),
                               ncol = length(sys_names),
                               sparse = TRUE)

        colnames(ssys) <- sys_names

        ssys_id <- mapply(match,
                          ssys_tmp,
                          MoreArgs = list(sys_names),
                          SIMPLIFY = FALSE)

        for (i in seq(along = ssys_id)) {
            if (length(ssys_id[[i]]) > 0) {
                ssys[i, ssys_id[[i]]] <- TRUE
            }
        }

    }
    else {
        ssys <- Matrix::Matrix(FALSE, nrow = nreact, ncol = 1, sparse = TRUE)
        colnames(ssys) <- NA
    }

    return(ssys)

}
