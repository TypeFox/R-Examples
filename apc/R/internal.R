#######################################################
#	apc package
#	Bent Nielsen, 1 Feb 2016, version 1.2
#	Internal functions
#######################################################
#	Copyright 2016 Bent Nielsen
#	Nuffield College, OX1 1NF, UK
#	bent.nielsen@nuffield.ox.ac.uk
#
#	This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#######################################################

apc.internal.function.date.2.character <- function(x,n.decimal=NULL)
{	#	apc.internal.function.date.2.character
	if(is.null(n.decimal))	return(as.character(x))
	if(n.decimal==1)		return(sprintf("%.1f",x))
	if(n.decimal==2)		return(sprintf("%.2f",x))
	if(n.decimal==3)		return(sprintf("%.3f",x))
	if(n.decimal==4)		return(sprintf("%.4f",x))
	if(n.decimal==5)		return(sprintf("%.5f",x))
	if(n.decimal==6)		return(sprintf("%.6f",x))
	if(n.decimal>6)			return(as.character(x))	
}	#	apc.internal.function.date.2.character
