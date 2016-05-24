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
"pHslope" <- 
function(Etris=-72.4,Eamp=4.9,S=35,T=25){
#	Etris=-99.4; Eamp=-18.9; S=35.1; T=1.6
	R  =  8.31447215 # J K–1 mol–1, Dickson et al. (2007)
	F  =  96485.339924 # C mol–1, Dickson et al. (2007)
	pHslope = 100*(((Eamp/1000-Etris/1000)/(tris(S=S,T=T)-amp(S=S,T=T)))-(R*(T + 273.15)*log(10)/F))/(R*(T + 273.15)*log(10)/F)
	attr(pHslope,"unit") <- "percent"
	return(pHslope)
}