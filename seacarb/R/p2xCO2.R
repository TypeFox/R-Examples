# Copyright (C) 2014 James Orr
# This file is part of seacarb.
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Convert pCO2 (uatm) to xCO2 (ppm) 

p2xCO2 <- function(S=35, T=25, Patm=1.0, pCO2=400){
  pH20 <- vapress(T=T, S=S, form="d2007")
  xCO2 <- pCO2 / (Patm - pH20) 
return(xCO2)
}
