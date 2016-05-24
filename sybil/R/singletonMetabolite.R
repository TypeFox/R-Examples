#  singletonMetabolite.R
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
# Function: .singletonMetabolite
#
#
#


.singletonMetabolite <- function(mat, exclM, exclR) {

    if ( (missing(exclM)) || (any(is.na(exclM))) ) {
        deCandM <- logical(nrow(mat))
    }
    else {
        deCandM <- exclM
    }

    if ( (missing(exclR)) || (any(is.na(exclR))) ) {
        deCandR <- logical(ncol(mat))
    }
    else {
        deCandR <- exclR
    }

    # generate smaller matrix without singleton metabolites
    indMatM <- c(1:nrow(mat))[!deCandM]
    indMatR <- c(1:ncol(mat))[!deCandR]

    tmp_mat <- mat[!deCandM, , drop = FALSE]
    tmp_mat <- tmp_mat[ , !deCandR, drop = FALSE]

    check  <- TRUE
    sreact <- logical(ncol(tmp_mat))
    smet   <- logical(nrow(tmp_mat))

    while(isTRUE(check)) {

        rs    <- rowSums(tmp_mat)

        # row indices of reactions used only once in S
        rrs <- which(rs == 1)
        smet[rrs] <- TRUE

        if (length(rrs) > 0) {
            # get reactions (columns) using singleton metabolites
            crs <- which(colSums(tmp_mat[rrs, , drop = FALSE]) != 0)
            tmp_mat[ , crs] <- FALSE
            sreact[crs] <- TRUE
        }
        else {
            check <- FALSE
        }
    }

    SMret <- logical(nrow(mat))
    SRret <- logical(ncol(mat))

    SMret[indMatM] <- smet
    SRret[indMatR] <- sreact

    return(list(smet = SMret, sreact = SRret))

}
