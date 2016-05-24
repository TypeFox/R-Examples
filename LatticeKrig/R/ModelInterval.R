
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

## LKrig model for 1-d data on an interval

LKrigSetupLattice.LKInterval <- function(object, x, verbose,
                                       NC, NC.buffer=5,  ...){
#object is usually of class LKinfo
  range.x<- rangeLocations<-  range(x)
  if( ncol( x) !=1) {
       stop( "x is not 1-d !")
   }
  grid.info<-list( xmin=  rangeLocations[1], xmax= rangeLocations[2],
   range= rbind(range.x))
  nlevel<- object$nlevel
  delta.level1<- ( grid.info$xmax - grid.info$xmin ) /( NC - 1 )
  mLevel<- rep( NA, nlevel)
  delta.save<- rep( NA, nlevel)
  grid.all.levels<- NULL
# begin multiresolution loop 
   for (j in 1:nlevel) {
        delta <- delta.level1/(2^(j - 1))
        delta.save[j] <- delta
        # the width in the spatial coordinates for NC.buffer grid points at this level.
        buffer.width <- NC.buffer * delta
        grid.list <- list(x = seq(grid.info$xmin - buffer.width, 
            grid.info$xmax + buffer.width, delta) )
        class( grid.list) <- "gridList"
        mLevel[j] <- length(grid.list$x)
        grid.all.levels <- c(grid.all.levels, list(grid.list))
    } 
# end multiresolution level loop
# create a useful index that indicates where each level starts in a
# stacked vector of the basis function coefficients.
    offset <- as.integer(c(0, cumsum(mLevel)))
    m<- sum(mLevel)
    mx<- cbind(mLevel)
    mLevelDomain <- (mLevel - 2 * NC)
    mxDomain<- cbind( mLevelDomain)
    out<-  list(  m = m, offset = offset, mLevel = mLevel,
                  delta = delta.save, rangeLocations = rangeLocations,
  # specific arguments for LKInterval              
                  mLevelDomain = mLevelDomain,
                  mx= mx, mxDomain= mxDomain,
                  NC=NC,
                  NC.buffer = NC.buffer,
                  grid = grid.all.levels,
                  grid.info=grid.info)
 return( out )
}

LKrigSAR.LKInterval<- function(object, Level, ... ){
   m<- object$latticeInfo$mLevel[Level] 
   a.wght<- (object$a.wght)[[Level]]
   if( length(a.wght) > 1) {
     stop("a.wght must be constant")
   }
   da<- c( m,m)
   ra<- c(rep( a.wght, m), rep( -1, (m-1)), rep( -1, (m-1)) )
   Bi <-  c( 1:m,2:m, 1:(m-1))
   Bj<- c( 1:m, 1:(m-1), 2:m)
  return(list(ind = cbind(Bi, Bj), ra = ra, da = da)) 
}  

LKrigLatticeCenters.LKInterval<- function(object, Level, ... ){
# return locations as a gridList object to be consistent with the
# Rectangle and Box geometries
   gridl<- object$latticeInfo$grid[[Level]]
# this would also work, however, if the grid locations were returned
#  explicitly as a one column matrix   
#    gridl<- matrix( unlist(object$latticeInfo$grid[[Level]]), ncol=1)
   return( gridl )
} 


  