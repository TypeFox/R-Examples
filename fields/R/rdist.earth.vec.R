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
rdist.earth.vec = function(x1, x2, miles=TRUE, R=NULL) {
  
  #set default radius
  if(is.null(R)) {
    if(miles)
      R = 3963.34
    else
      R = 6378.388
  }
  
  #convert lon/lat to radians
  x1 = x1 * (pi/180)
  x2 = x2 * (pi/180)
  
  #calculate distances using Haversine method
  lonDist2 = (x2[,1] - x1[,1]) * (1/2)
  latDist2 = (x2[,2] - x1[,2]) * (1/2)
  a = sin(latDist2) * sin(latDist2) + cos(x1[, 2]) * cos(x2[, 2]) * sin(lonDist2) * sin(lonDist2)
  return(2 * atan2(sqrt(a), sqrt(1 - a)) * R)
}
