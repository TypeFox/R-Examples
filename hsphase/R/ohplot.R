# Copyright (C) 2014 Mohammad H. Ferdosi and Cedric Gondro  
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

ohplot <- function(oh, genotype, pedigree, check = FALSE)
{
    if (check)
    {
		if(!is.matrix(oh))
		{
			stop("The oh must be matrix")
		}
        if (any(duplicated(rownames(genotype))))
        {
            stop("There are duplicated ID in the input file")
        }
        if (is.null(genotype)) 
            stop("Invalid input!")
        if (!is.matrix(genotype)) 
            stop("GenotypeMatrix should be a MATRIX")
        if (length(genotype[genotype != 9 & genotype != 0 & genotype != 2 & genotype != 1]) > 0) 
            stop("GenotypeMatrix must contain only 0, 1, 2 or 9")
		if (ncol(pedigree) < 2) 
			stop("Pedigree must have at least 2 columns")
		if (any(duplicated(pedigree[,1])))
		{
			stop("There are duplicated ID in the input file")
		}
		if(nrow(genotype)!=nrow(oh))
		{
			stop("The number of rows in the oh and genotype must be equal") 
		}
		if(ncol(oh)!=nrow(oh))
		{
			stop("The oh matrix mast have the same number of rows and columns")
		}
    }
    
	genotype[is.na(genotype)] <- 9	
	sorted <- sort(oh[upper.tri(oh)])
	maxOH <- max(length(sorted))
	difs <- sorted[2:length(sorted)] - sorted[1:(length(sorted) - 1)]
	index <- which(difs == max(difs))
	bestsep <- sorted[index + 1] - sorted[index]
	sepval <- round(bestsep/max(sorted), 3)
	cutoff <- round((sorted[index + 1] + sorted[index])/2)
	
	
	if (missing(pedigree))
	{
		plot(sorted, main = paste("Separation value:", sepval, "\ncut off (number of SNP):", cutoff), pch = 20, col = "blue", xlab = "", ylab = "opposing homozygotes", 
				xaxt = "n")
		p <- apply(genotype, 2, .maf)
		p <- p[!is.na(p)]
		q <- 1 - p
		halfsibs <- sum(p^2 * q^2)
		unrelated <- sum(2 * (p^2 * q^2))
		fullsib <- sum(0.5 * (p^2 * q^2))
		maxsnpnooh <- (sum(p^2 * q^2) + sum(2 * (p^2 * q^2)))/2
		maxsnpnooh <- maxsnpnooh - (0.1 * maxsnpnooh)
		
		abline(h = c(cutoff, halfsibs, unrelated, fullsib, maxsnpnooh), col = 2:6, lty = 2:6)
		legend("bottomright", legend = c("cut off", "average half-sib", "average unrelated", "average full-sib", "90% of average"), col = c(2, 3, 4, 
						5, 6), lty = 2:6, cex = 0.7)
	} else
	{
		colorMatrix <- matrix(2, nrow = nrow(oh), ncol = ncol(oh))
		rownames(colorMatrix) <- rownames(oh)
		colnames(colorMatrix) <- colnames(oh)
		pedigree <- as.matrix(pedigree)
		ped <- pedigree[pedigree[, 1] %in% rownames(colorMatrix), ]
		ped <- pedigree[pedigree[, 2] %in% rownames(colorMatrix), ]
		if(nrow(ped)>4)
		{
			colorMatrix[ped] <- 1
			colorMatrix[ped[, c(2, 1)]] <- 1
		}
		orderCol <- order(as.vector(oh[upper.tri(oh)]))
		halfsibsFamily <- tapply(pedigree[, 1], pedigree[, 2], function(x) x)
		halfsibsFamily <- lapply(halfsibsFamily, function(x) expand.grid(x, x))
		halfsibsFamily <- as.matrix(do.call(rbind, halfsibsFamily))
		if(nrow(halfsibsFamily)>1)
		{
			colorMatrix[halfsibsFamily[, c(2, 1)]] <- 4
			colorMatrix[halfsibsFamily] <- 4
		}
		
		plot(oh[upper.tri(oh)][orderCol], main = paste("separation value:", sepval, "\ncut off (number of SNP):", cutoff), pch = 20, col = colorMatrix[upper.tri(colorMatrix)][orderCol], 
				xlab = "", ylab = "opposing homozygotes", xaxt = "n")
		p <- apply(genotype, 2, .maf)
		p <- p[!is.na(p)]
		q <- 1 - p
		halfsibs <- sum(p^2 * q^2)
		unrelated <- sum(2 * (p^2 * q^2))
		fullsib <- sum(0.5 * (p^2 * q^2))
		maxsnpnooh <- (sum(p^2 * q^2) + sum(2 * (p^2 * q^2)))/2
		maxsnpnooh <- maxsnpnooh - (0.1 * maxsnpnooh)
		
		abline(h = c(cutoff, halfsibs, unrelated, fullsib, maxsnpnooh), col = 2:6, lty = 2:6)
		legend("bottomright", legend = c("cut off", "average half-sib", "average unrelated", "average full-sib", "90% of average"), col = c(2, 3, 4, 
						5, 6), lty = 2:6, cex = 0.7)
		legend("topleft", legend = c("parent-offspring", "unrelated", "half-sib"), fill = c(1, 2, 4), cex = 0.7)
	}
} 
