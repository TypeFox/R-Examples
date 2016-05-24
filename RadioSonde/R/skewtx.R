"skewtx" <-
function(temp, ycoord)
{
#
# Copyright 2001,2002 Tim Hoar
#
# This file is part of the RadioSonde library for R and related languages.
#
# RadioSonde is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# RadioSonde is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RadioSonde; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
	#
	# Determine x-coordinate on skew-T, log p diagram given
	#           temperature (C)
	# and y-coordinate from FUNCTION SKEWTY.  X-origin at T=0c.
	#
	0.54000000000000004 * temp + 0.90691999999999995 * ycoord
}
