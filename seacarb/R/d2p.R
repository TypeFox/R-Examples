# Copyright (C) 2010 Jean-Pierre Gattuso and Héloïse Lavigne
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


#Reference:  
#    Saunders, Peter M. (1981) Practical conversion of pressure
#     to depth, J. Phys. Oceanogr., 11, 573-574, (1981)
#
#     From FORTRAN code of

#     R. Millard
#     March 9, 1983
#     check value: p80=7500.004 dbars;for lat=30 deg., depth=7321.45 meters

#     Modified (slight format changes + added ref. details):
#     J. Orr, MEL-IAEA, 16 April 2009



"d2p" <- function (depth, lat=40)

{
## Checking dimension of input data
ndepth <- length(depth)
nlat <- length(lat)
nmax <- max(ndepth, nlat)
if(ndepth < nmax){depth <- rep(depth[1], nmax)}
if(nlat < nmax){lat <- rep(lat[1], nmax)}

pi <- 3.141592654
pressure <- rep(NA,nmax)

for(i in (1:nmax)){

## Checking for lat positive
lat[i] <- abs(lat[i])

plat <- abs(lat[i]*pi)/180 
d <- sin(plat)
c1 <- 5.92e-3 + d^2 * 5.25e-3

pressure[i] <- ((1-c1)-sqrt(((1-c1)^2)-(8.84e-6*depth[i]))) / 4.42e-6

}

return(pressure)

}

