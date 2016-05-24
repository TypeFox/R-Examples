#  deadEndMetabolite.R
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
# Function: .deadEndMetabolite
#
#
# The function .deadEndMetabolite() is inspired by the functions
# detectDeadEnds() and removeDeadEnds contained in the COBRA Toolbox.
# The algorithm is basically the same.


.deadEndMetabolite <- function(mat,
                               lb,
                               exclM,
                               exclR,
                               tol = SYBIL_SETTINGS("TOLERANCE")) {

    # candM and candR can be used for singleton metabolites
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

    check   <- TRUE
    DEreact <- logical(ncol(tmp_mat))
    DEmet   <- logical(nrow(tmp_mat))


    # check is set to FALSE, if no more dead end metabolites can be found
    while (isTRUE(check)) {

        temp_dem <- rep(FALSE, nrow(tmp_mat))

        # metabolites
        i <- 1
        while(i <= nrow(tmp_mat)) {

            # isDE is FALSE, if metabolite i is not dead end
            isDE  <- NA

            # reactions (non-zero elements in row i)
            nz <- which(abs(tmp_mat[i, ]) > tol)

            j <- 1
            while(j <= length(nz)) {

                if ( (!is.na(isDE)) && (isTRUE(isDE)) ) {
                    break
                }

                if (is.na(isDE)) {

                    # assume, metabolite i is dead end for reaction j
                    isDE <- TRUE

                    # check, if metabolite i can be excluded from the
                    # dead end list: test all reactions k > j
                    k <- j+1
                    while (k <= length(nz)) {
                        if ( ( sign(tmp_mat[i, nz[k]]) != sign(tmp_mat[i, nz[j]]) ) ||
                             (lb[nz[k]] < 0)        ||
                             (lb[nz[j]] < 0) ) {
                            #print(paste(i,":",k))
                            isDE <- FALSE
                            break
                        }
                        else {
                            k <- k+1
                        }
                    }
                }

                # if isDE stays TRUE, metabolite i is dead end
                if (isTRUE(isDE)) {
                    temp_dem[i] <- TRUE
                    break
                }
                else {
                    j <- j + 1
                }
            }

            i <- i+1

        }

        if (sum(temp_dem) == 0) {
            check <- FALSE
        }
        else {
            DEmet[temp_dem] <- TRUE
            matb <- abs(tmp_mat) > tol
            crs  <- which(colSums(matb[temp_dem, , drop = FALSE]) != 0)
            tmp_mat[ , crs] <- 0
            DEreact[crs] <- TRUE
        }

    }


    # put DEmet and DEreact to their original size
    # (if arguments exclR or exclM were used)
    DEMret <- logical(nrow(mat))
    DERret <- logical(ncol(mat))

    DEMret[indMatM] <- DEmet
    DERret[indMatR] <- DEreact

    return(list(dem = DEMret, der = DERret))

}
