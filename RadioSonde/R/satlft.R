"satlft" <-
function(thw, p)
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

	cta <- 273.14999999999998
	akap <- 0.28541
	pwrp <- (p/1000.)^akap
	tone <- (thw + cta) * pwrp - cta
	eone <- wobf(tone) - wobf(thw)
	rate <- 1.
	dlt <- 1.
	while(abs(dlt) > 0.10000000000000001) {
		ttwo <- tone - eone * rate
		pt <- (ttwo + cta)/pwrp - cta
		etwo <- pt + wobf(ttwo) - wobf(pt) - thw
		dlt <- etwo * rate
		rate <- (ttwo - tone)/(etwo - eone)
		tone <- ttwo
		eone <- etwo
	}
	ttwo - dlt
}
