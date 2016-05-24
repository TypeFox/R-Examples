# Copyright (C) 2009 Jean-Pierre Gattuso
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
"amp" <- 
function(S=35,T=25){
	amp = (111.35+5.44875*S)*1/(T+273.15)+41.6775-0.015683*S-6.20815*log(T+273.15)-log10(1-0.00106*S)
	attr(amp,"unit") <- "mol/kg"
	return(amp)
}