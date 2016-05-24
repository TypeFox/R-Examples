# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKrig.precision <- function(LKinfo, return.B = FALSE,
                                   verbose=FALSE) { 
    L <- LKinfo$nlevel
    offset <- LKinfo$latticeInfo$offset
    # some checks on arguments
    LKinfoCheck(LKinfo)
    # ind holds non-zero indices and ra holds the values
    ind <- NULL
    ra <- NULL
    da <- rep(0, 2)
        # loop over levels
        for (j in 1:L) {
            # evaluate the SAR matrix at level j.
            tempB<- LKrigSAR( LKinfo, Level=j)
            if( verbose){
            	cat("dim indices in spind of B:",dim( tempB$ind) , fill=TRUE)            	
            }
            # multiply this block by 1/ sqrt(diag( alpha[[j]]))
            alpha.level <- (LKinfo$alpha)[[j]]
            if( verbose){
                cat("length alpha parameter", length( alpha.level), fill=TRUE)
            }
            if( any( is.na(c(alpha.level)) ) ){
                	stop("NAs in alpha list")}             	
            if (length(alpha.level) == 1) {
                tempra <- 1/sqrt(alpha.level[1]) * tempB$ra
            }
            else {
                rowindices <- tempB$ind[, 1]
                tempra <- 1/sqrt(alpha.level[rowindices]) * tempB$ra
            }
            # accumulate the new block
            # for the indices that are not zero
            ra <- c(ra, tempra)
            ind <- rbind(ind, tempB$ind + offset[j])
            # increment the dimensions
            da[1] <- da[1] + tempB$da[1]
            da[2] <- da[2] + tempB$da[2]
        }
        # dimensions of the full matrix
        # should be da after loop
        # check this against indices in LKinfo
        #
        if ((da[1] != offset[L + 1]) | (da[2] != offset[L + 
            1])) {
            stop("Mismatch of dimension with size in LKinfo")
        }
        # convert to spind format:
        # tempBtest <- list(ind = ind, ra = ra, da = da) 
        # tempB<- spind2spam( tempBtest)
        if( verbose){
        	cat("dim of ind (fullB):", dim( ind), fill=TRUE)
        }
        tempB <- spam( list( ind=ind, ra), nrow=da[1], ncol=da[2])
         if( verbose){
        	cat("dim after spind to spam in precision:", dim( tempB), fill=TRUE)
        }
    if (return.B) {
        return(tempB)
    }
    else {
        # find precision matrix Q = t(B)%*%B and return
        return(t(tempB) %*% (tempB))
    }
}

