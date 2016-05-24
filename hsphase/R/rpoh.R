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

rpoh <- function(genotypeMatrix, oh, forwardVectorSize = 30,excludeFP = TRUE,nsap = 3, maxRec = 15, intercept = 26.3415, coefficient = 77.3171, snpnooh, method, maxsnpnooh)
{
	
	if(missing(method))
		stop("please set the method")
	METHODS <- c("recombinations", "simple","calus","manual" )
	method <- pmatch(method, METHODS)
	
	if (is.na(method)) 
		stop("invalid pedigree reconstruction method")
	if (method == -1) 
		stop("ambiguous pedigree reconstruction method")
	if(method == 1)
	{
		if(forwardVectorSize == 30 && excludeFP == TRUE && nsap == 3 && maxRec == 15)
		{
			print("Default values will be used for identification of recombinations")
		}
		result <- .rpohhsphase(genotypeMatrix = genotypeMatrix,oh = oh,forwardVectorSize = forwardVectorSize,excludeFP = excludeFP,nsap = nsap, maxRec = maxRec)	
	}
	if(method == 2)
	{
		if(intercept == 26.3415 && coefficient == 77.3171)
		{
			print("Default values will be use for caculation of maximum possible recombinations")
			print(paste(snpnooh * 1000, "SNPs were used to create the oh matrix"))
		}
		result <- .prSimple(oh, snpnooh, intercept = intercept, coefficient = coefficient) 	
	}
	if(method == 3)
	{
		result <- .prCalus(oh, genotypeMatrix)
	}
	if(method == 4)
	{
		result <- .prManual(oh, maxsnpnooh)
	}
	result[,2] <- as.factor(result[,2])
	levels(result[,2]) <- 1:length(levels(result[,2]))
	result
}

