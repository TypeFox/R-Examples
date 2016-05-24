# Copyright (C) 2014 Mohammad H. Ferdosi
#
# HSPhase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# HSPhase program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http:#www.gnu.org/licenses/>.
ssp <- function(blockMatrix, genotypeMatrix)
{
    if (nrow(blockMatrix) != nrow(genotypeMatrix)) 
        stop("blockMatrix and genotypeMatrix must have the same number of rows")
    if (!is.matrix(blockMatrix)) 
        stop("blockMatrix should be a MATRIX")
    if (!is.matrix(genotypeMatrix)) 
        stop("genotypeMatrix should be a MATRIX")
    if (length(genotypeMatrix[genotypeMatrix != 0 & genotypeMatrix != 2 & genotypeMatrix != 1 & genotypeMatrix != 
        9]) > 0) 
        stop("genotypeMatrix must contain only 0,1 and 2")
    if (length(blockMatrix[blockMatrix != 0 & blockMatrix != 1 & blockMatrix != 2]) > 0) 
        stop("blockMatrix must contain only 0, 1 and 2")
    expandMat <- as.numeric(genotypeMatrix)
    n <- ncol(genotypeMatrix)
    fMat <- t(matrix(as.double(rep(0, n * 2)), ncol = ncol(genotypeMatrix)))
    siregenotype <- apply(genotypeMatrix, 2, function(x)
    {
        if (any(x == 0) && any(x == 2))
        {
            x <- 1
        }
        else
        {
            x <- 0
        }
    })
    result <- .C("ssp", block = as.integer(blockMatrix), genotype = as.integer(genotypeMatrix), siregenotype = as.integer(siregenotype), 
        nrow = as.integer(nrow(genotypeMatrix)), ncol = as.integer(ncol(genotypeMatrix)), result = fMat)$result
    result <- t((result))
	if(!is.null(colnames(genotypeMatrix)))
	    colnames(result) <- colnames(genotypeMatrix)
	result
} 
