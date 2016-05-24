# $Id: kinresplot.R 106 2011-05-12 10:55:37Z jranke $

# Copyright (C) 2008-2010 Johannes Ranke
# Contact: mkin-devel@lists.berlios.de

# This file is part of the R package kinfit

# kinfit is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

kinresplot <- function(kinobject, kinmodel,
	xlab = "Time [days]", ylab = "Residual [% of applied radioactivity]",
	maxabs = "auto")
{
	m <- kinobject$fits[[kinmodel]]
	t <- m$model$t
	residuals <- residuals(m)
	if (maxabs == "auto") maxabs = max(abs(residuals), na.rm = TRUE)
	plot(t, residuals,
		xlab = xlab,
		ylab = ylab,
		ylim = c( -1.2 * maxabs, 1.2 * maxabs))
  abline(h=0, lty=2)
	title(paste("Residuals of", kinmodel, "fit"))
}
