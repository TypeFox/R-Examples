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

print.LKinfo <- function(x, ...) {
    LKinfo <- x
    L <- LKinfo$nlevel
    cat("Classes for this object are: " , class( LKinfo), fill=TRUE)
    cat("The second class usually will indicate the geometry
     e.g.  2-d rectangle is  LKRectangle", fill=TRUE)
    cat(" ", fill = TRUE) 
    cat("Ranges of locations in raw scale:", fill=TRUE)
    print(  LKinfo$latticeInfo$rangeLocations)
    if( !is.null(LKinfo$basisInfo$V)){
    	cat("(inverse) linear transformation for lattice nodes:",fill=TRUE)
    	print(LKinfo$basisInfo$V )
    	cat("transformed ranges:",fill=TRUE)
    	print( LKinfo$latticeInfo$grid.info$range)
    }
    cat(" ", fill = TRUE)
    cat("Number of levels:", L, fill = TRUE)
    cat("delta scalings:", x$latticeInfo$delta, fill = TRUE)
    cat("with an overlap parameter of ", LKinfo$basisInfo$overlap, fill=TRUE)
    cat("alpha: ", unlist(x$alpha), fill = TRUE)
    if (!is.null(x$nu)) {
        cat("based on smoothness nu = ", x$nu, fill = TRUE)
    }
    cat("a.wght: ", unlist(x$a.wght), fill = TRUE)
       cat(" ", fill = TRUE)
  # Details on basis functions at each level
      bType <- LKinfo$basisInfo$BasisType
      cat( "Basis  type:",
           LKinfo$basisInfo$BasisType, 
           "using ",
           LKinfo$basisInfo$BasisFunction,
           " and", LKinfo$distance.type, " distance.", fill=TRUE)
               if( LKinfo$normalize){
  cat("Basis functions will be normalized", fill=TRUE)
        }
  cat(" ", fill = TRUE)      
  cat("Total number of basis functions ",  LKinfo$latticeInfo$m, fill=TRUE)  
      temp<- cbind( 1:L, LKinfo$latticeInfo$mLevel)
      dimnames( temp)<- list( rep( "", L), c("Level" ,   "Basis size"))
     if( !is.null(LKinfo$latticeInfo$mx)){
    	 temp<- cbind( temp, LKinfo$latticeInfo$mx)
 	 }
      print( temp)      	
        cat(" ", fill = TRUE)  
 cat("Lambda value: ", LKinfo$lambda, fill=TRUE)         
}



