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

aio <- function(genotypeMatrix, bmh_forwardVectorSize = 30, bmh_excludeFP = TRUE, bmh_nsap = 3, output = "phase")
{
    block <- bmh(genotypeMatrix, forwardVectorSize = bmh_forwardVectorSize, excludeFP = bmh_excludeFP, nsap = bmh_nsap)
    id <- rownames(genotypeMatrix)
    id1 <- paste(id, "p", sep = "")
    id2 <- paste(id, "m", sep = "")
    id <- character(length(id) * 2)
    id[seq(1, nrow(genotypeMatrix) * 2, 2)] <- id1
    id[seq(2, nrow(genotypeMatrix) * 2, 2)] <- id2
    sireHaplotype <- ssp(block, genotypeMatrix)
    singleStrand <- phf(genotypeMatrix, block, sireHaplotype)
    phaseResult <- matrix(rep(9, ncol(genotypeMatrix) * nrow(genotypeMatrix) * 2), ncol = ncol(genotypeMatrix))
    st1 <- seq(from = 1, to = nrow(genotypeMatrix) * 2, by = 2)
    st2 <- seq(from = 2, to = nrow(genotypeMatrix) * 2, by = 2)
    phaseResult[st1, ] <- singleStrand
    phaseResult[st2, ] <- genotypeMatrix - singleStrand
    phaseResult[phaseResult != 0 & phaseResult != 1] <- 9
    if (!is.null(rownames(genotypeMatrix)))
    {
        rownames(phaseResult) <- id
    }
    colnames(phaseResult) <- colnames(genotypeMatrix)
    if (output != "phase") 
        list(phasedHalfsibs = phaseResult, sireHaplotype = sireHaplotype, blockStructure = block)
    else phaseResult
} 
