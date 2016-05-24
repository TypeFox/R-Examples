# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
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
compactToMat = function(compactMat, diagVal=0, lower.tri=FALSE, upper.tri=TRUE) {
  #compactMat: a symmetric matrix stored as a vector containing elements for the upper triangle
  #portion of the true matrix
  #diagVal: a number to put in the diagonal entries of the output matrix
  #lower.tri: if TRUE, fills in lower tringular portion of the matrix
  #upper.tri: if TRUE, fills in upper tringular portion of the matrix
  
  if(class(compactMat) == 'dist') {
    n <- attr(compactMat, "Size")
  } else { # (n^2 - n)/2 = length(compactMat)
    stop("input matrix is not compact or is not of class \"dist\"")
    
    #or if class is not dist but input matrix is still compact, use:
    #n = (1 + sqrt(1 + 8*length(compactMat)))/2
  }
  
  return(.Call("compactToMatC", as.double(compactMat),
               as.integer(length(compactMat)), 
               as.integer(n), as.double(diagVal), 
               as.integer(lower.tri), as.integer(upper.tri)))
}
