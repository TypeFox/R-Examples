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

impute <- function(halfsib_genotype_ld, sire_hd,bmh_forwardVectorSize = 30, bmh_excludeFP = TRUE, bmh_nsap = 3)
{
	if (length(colnames(halfsib_genotype_ld)) == 0) 
		stop("The halfsib_genotype_ld must have a colnames")
	
	if (length(colnames(sire_hd)) == 0) 
		stop("The sire_hd must have a colnames")
	
	halfsib_block_ld <- bmh(halfsib_genotype_ld,forwardVectorSize = bmh_forwardVectorSize, excludeFP = bmh_excludeFP, nsap = bmh_nsap)
	sire_ld <- ssp(halfsib_block_ld, halfsib_genotype_ld)
	halfsib_paternal_ld <- phf(halfsib_genotype_ld, halfsib_block_ld, sire_ld)
	
	snpld <- colnames(halfsib_paternal_ld)
	snphd <- colnames(sire_hd)
	ldindex <- match(snpld, snphd)
	nald <- which(is.na(ldindex))
	if (length(nald) > 0)
	{
		ldindex <- ldindex[-nald]
		halfsib_paternal_ld <- halfsib_paternal_ld[, -nald]
		halfsib_block_ld <- halfsib_block_ld[, -nald]
		sire_ld <- sire_ld[, -nald]
	}
	
	sire_hdformLD <- sire_hd[,ldindex]
	
	sire_hdformLD[sire_hdformLD==9] <- NA
	sire_ld[sire_ld==9] <- NA
	a <- summary(lm(sire_ld[1, ] ~ sire_hdformLD[1, ],na.action = "na.omit"))$r.squared
	b <- summary(lm(sire_ld[1, ] ~ sire_hdformLD[2, ],na.action = "na.omit"))$r.squared
	
	if (a < b)
	{
		sire_hd <- sire_hd[c(2, 1), ]
	}
	
	result <- matrix(0, ncol = ncol(sire_hd), nrow = nrow(halfsib_paternal_ld))
	result[, ldindex] <- halfsib_block_ld
	result[result == 1] <- 3
	result[result == 2] <- 4
	result <- .hblock(result)
	result <- .phfnoGenotype(result, sire_hd)
	rownames(result) <- rownames(halfsib_genotype_ld)
	colnames(result) <- colnames(sire_hd)
	result
} 
