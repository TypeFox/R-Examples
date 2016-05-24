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
"tris" <- 
function(S=35,T=25){
	tris = ((11911.08-18.2499*S-0.039336*S^2)*1/(T+273.15))-366.27059+0.53993607*S+0.00016329*S^2+((64.52243-0.084041*S)*log(T+273.15))-(0.11149858*(T+273.15))
	attr(tris,"unit") <- "mol/kg"
	return(tris)
}