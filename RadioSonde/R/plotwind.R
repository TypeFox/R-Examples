"plotwind" <-
function(dataframe, size = 5., ylim = c(1050, 100), legend = FALSE)
{
#
# Copyright 2001,2002 Eric Gilleland, and Doug Nychka
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
	y <- dataframe$press
	y <- y[!is.na(y)]
	y <- skewty(y)
	ylim <- skewty(ylim)
	n <- length(y)
	u <- dataframe$uwind
	v <- dataframe$vwind
	theta <- dataframe$dir
	speed <- dataframe$wspd
	speed <- round(speed, digits = 0.)
	speed <- as.integer(speed)
	plot(c(-1.5, 1.5), ylim, axes = FALSE, type = "n", xlab = "", ylab = "")
	if(legend) {
		mtext("full barb = 10 m/s", side = 1, line = 1)
	}
	abline(v=0)
	for(k in 1:length(y)) {
		if(y[k] > ylim[2])
			break
		if(!is.na(speed[k])) {
			if(speed[k] != 999.) {
	station.symbol(0., y[k], speed = speed[k], direction = 
					theta[k], circle=FALSE, cex = size)
			}
		}
	}
}
