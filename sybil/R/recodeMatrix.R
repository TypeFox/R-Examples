#  recodeMatrix.R
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
# Function: .recodeMatrix
#
# arguments:
#
# mat:   the matrix to recode (class Matrix)
# signs: character vector of length three:
#        mat[i, j]  > tol    -->   signs[3]
#        mat[i, j]  < tol    -->   signs[1]
#        mat[i, j] == tol    -->   signs[2]
# tol:   tolerance value


.recodeMatrix <- function(mat,
                          signs = c("-", " ", "+"),
                          tol = SYBIL_SETTINGS("TOLERANCE")) {

    stopifnot(is(mat, "Matrix"), length(signs) == 3)

    mat <- apply(mat, 2,
                 function(x) {
                     ifelse(x > tol, signs[3], ifelse(abs(x) > tol,
                                                      signs[1], signs[2]))
                 })

    return(mat)

}
