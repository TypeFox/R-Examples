# Copyright (C) 2013 Mohammad H. Ferdosi
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


.hblock <- function(bmhResult, MaxBlock = 400)
{
    if (!is.matrix(bmhResult)) 
        stop("GenotypeMatrix should be a MATRIX")
    expandMat <- as.double(t(bmhResult))
    n <- ncol(bmhResult) * nrow(bmhResult)
    fMat <- matrix(as.integer(rep(0, n)), nrow = ncol(bmhResult))
    result <- .C("hblock", expandMat = as.integer(expandMat), nrow = as.integer(ncol(bmhResult)), ncol = as.integer(nrow(bmhResult)), 
        result = fMat, MB = as.integer(MaxBlock))$result
    colnames(result) <- rownames(bmhResult)

    if (!is.null(colnames(bmhResult))) 
        rownames(result) <- colnames(bmhResult)
    t(result)
}



hbp <- function(PhasedGenotypeMatrix, PhasedSireGenotype, strand = "auto")
{
	if (!is.matrix(PhasedGenotypeMatrix)) 
		stop("PhasedGenotypeMatrix should be a MATRIX")
	if (length(PhasedGenotypeMatrix[PhasedGenotypeMatrix != 0 & PhasedGenotypeMatrix != 1 & 
							PhasedGenotypeMatrix != 9]) > 0) 
		stop("PhasedGenotypeMatrix must contain only 0 and 1 or 9 for missing SNP")
	if (!is.matrix(PhasedSireGenotype)) 
		stop("PhasedSireGenotype should be a MATRIX")
	if (length(PhasedSireGenotype[PhasedSireGenotype != 0 & PhasedSireGenotype != 1 & PhasedSireGenotype != 
							9]) > 0) 
		stop("PhasedSireGenotype must contain only 0 and 1 or 9 for missing SNP")
	if (ncol(PhasedGenotypeMatrix) != ncol(PhasedSireGenotype)) 
		stop("Number of markers in sire and half-sib family must be the same")
	if (nrow(PhasedSireGenotype) != 2) 
		stop("PhasedSireGenotype must have 2 rows")
	
	METHODS <- c("auto", 1,2)
	method <- pmatch(strand, METHODS)
	if (is.na(method)) 
		stop("invalid strand")
	if (method == -1) 
		stop("ambiguous pedigree reconstruction method")
	if (method == 1) 
	{
		str = 0
	}
	else if(method == 2)
	{
		str = 1		
	}
	else if(method == 3)
	{
		str = 2
	}
	expandMat <- as.numeric(t(PhasedGenotypeMatrix))
	n <- ncol(PhasedGenotypeMatrix) * nrow(PhasedGenotypeMatrix)
	fMat <- matrix(as.integer(rep(0, n/2)), nrow = ncol(PhasedGenotypeMatrix))
	result <- .C("hbphased", expandMat = as.integer(expandMat), nrow = as.integer(nrow(PhasedGenotypeMatrix)), 
			ncol = as.integer(ncol(PhasedGenotypeMatrix)), result = fMat, siregenotype = as.integer(t(PhasedSireGenotype)), strand = as.integer(str))$result
	
	if (!is.null(rownames(PhasedGenotypeMatrix))) 
		colnames(result) <- rownames(PhasedGenotypeMatrix)[seq(from = 1, by = 2, to = nrow(PhasedGenotypeMatrix))]
	t(result)
}
phf <- function(GenotypeMatrix, blockMatrix, sirePhasedMatrix)
{
	if (is.null(GenotypeMatrix) | is.null(blockMatrix) | is.null(sirePhasedMatrix)) 
		stop("Invalid input!")
	if ((!is.matrix(GenotypeMatrix)) | (!is.matrix(blockMatrix)) | (!is.matrix(sirePhasedMatrix))) 
		stop("All inputs should be a MATRIX")
	if (length(GenotypeMatrix[GenotypeMatrix != 0 & GenotypeMatrix != 2 & GenotypeMatrix != 
							1 & GenotypeMatrix != 9]) > 0) 
		stop("GenotypeMatrix must contain only 0,1 and 2 or 9 for missing SNPs")
	if (length(blockMatrix[blockMatrix != 0 & blockMatrix != 2 & blockMatrix != 1]) > 0) 
		stop("blockMatrix must contain only 0,1 and 2")
	if (length(sirePhasedMatrix[sirePhasedMatrix != 0 & sirePhasedMatrix != 1 & sirePhasedMatrix != 
							9]) > 0) 
		stop("SireMatrix must contain only 0,1 and 9")
	expandMat <- as.numeric(GenotypeMatrix)
	n <- nrow(GenotypeMatrix) * ncol(GenotypeMatrix)
	fMat <- matrix(as.integer(rep(0, n)), ncol = ncol(GenotypeMatrix))
	result <- .C("phase", genotype = as.integer(GenotypeMatrix), nrow = as.integer(nrow(GenotypeMatrix)), 
			ncol = as.integer(ncol(GenotypeMatrix)), block = as.integer(blockMatrix), sirePhasedMatrix = as.integer(t(sirePhasedMatrix)), 
			result = fMat)$result
   colnames(result) <- colnames(GenotypeMatrix)
   result
}
.phfnoGenotype <- function(blockMatrix, sirePhasedMatrix)
{
		
	if (length(blockMatrix[blockMatrix != 0 & blockMatrix != 2 & blockMatrix != 1]) > 0) 
		stop("blockMatrix must contain only 0,1 and 2")
	if (length(sirePhasedMatrix[sirePhasedMatrix != 0 & sirePhasedMatrix != 1 & sirePhasedMatrix != 
							9]) > 0) 
		stop("SireMatrix must contain only 0,1 and 9")
	expandMat <- as.numeric(blockMatrix)
	n <- nrow(blockMatrix) * ncol(blockMatrix)
	fMat <- matrix(as.integer(rep(0, n)), ncol = ncol(blockMatrix))
	.C("phaseNogenotype", nrow = as.integer(nrow(blockMatrix)), 
			ncol = as.integer(ncol(blockMatrix)), block = as.integer(blockMatrix), sirePhasedMatrix = as.integer(t(sirePhasedMatrix)), 
			result = fMat)$result
}
pm <- function(blockMatrix, method = "constant")
{
	## if(missing(method))
	##     stop("please set the method")
	METHODS <- c("constant", "relative")
	method <- pmatch(method, METHODS)
	
	if (is.na(method)) 
		stop("invalid method")
	if (method == -1) 
		stop("ambiguous  method")
	
	if (!is.matrix(blockMatrix)) 
		stop("blockMatrix should be a MATRIX")
	expandMat <- as.double(t(blockMatrix))
	n <- ncol(blockMatrix) * nrow(blockMatrix)
	fMat <- matrix(as.double(rep(0, n)), nrow = ncol(blockMatrix))
	res <- .C("pm", expandMat = as.integer(expandMat), nrow = as.integer(ncol(blockMatrix)), 
			ncol = as.integer(nrow(blockMatrix)), method = method,result = fMat)$result
	result <- t(res[-nrow(res), ])
	rownames(result) <- rownames(blockMatrix)
	colnames(result) <- colnames(blockMatrix)[-1]
	result
}


recombinations <- function(blockMatrix)
{
    if (!is.matrix(blockMatrix)) 
        stop("The inputs must be MATRIX")
    if (length(blockMatrix[blockMatrix != 0 & blockMatrix != 1 & blockMatrix != 2]) > 0) 
        stop("Inputs must contain only 0 and 1 or 2")
    mat <- as.numeric(t(blockMatrix))
    nSwitch <- integer(nrow(blockMatrix))
    .C("recombinations", mat = as.integer(mat), nrow = as.integer(nrow(blockMatrix)), ncol = as.integer(ncol(blockMatrix)), 
        result = nSwitch)$result
} 


