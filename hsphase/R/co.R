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

co <- function(genotypeMatrix)
{
	if (is.null(genotypeMatrix)) 
		stop("Invalid input!")
	if (!is.matrix(genotypeMatrix)) 
		stop("genotypeMatrix should be a MATRIX")
	if (length(genotypeMatrix[genotypeMatrix != 9 & genotypeMatrix != 0 & genotypeMatrix != 2 & genotypeMatrix != 1]) > 0) 
		stop("genotypeMatrix must contain only 0, 1, 2 or 9")
	
	
	zeroFreq <- colMeans(genotypeMatrix == 0) * nrow(genotypeMatrix)
	oneFreq <- colMeans(genotypeMatrix == 1) * nrow(genotypeMatrix)
	twoFreq <- colMeans(genotypeMatrix == 2) * nrow(genotypeMatrix)
	
	hetsite <- apply(genotypeMatrix, 2, function(x)
			{
				if (any(x == 0) && any(x == 2))# && any(x==1))
				{
					x <- 1
				} else
				{
					x <- 0
				}
			})
	
	genotypeMatrix[,which(hetsite==0)] <- 9;
	genotypeMatrix[genotypeMatrix==1] <- 9;
	result <- .Call("co", genotypeMatrix, hetsite, PACKAGE = "hsphase")
	result[,-ncol(result)]
} 
