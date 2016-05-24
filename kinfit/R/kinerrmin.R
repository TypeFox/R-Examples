# $Id: kinerrmin.R 59 2010-07-28 12:29:15Z jranke $

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

kinerrmin <- function(kinfits, kinmodel = "SFO", alpha = 0.05)
{
	m = kinfits[[kinmodel]]

	kindata <- data.frame(t = kinfits[[kinmodel]]$model$t, 
		parent = kinfits[[kinmodel]]$model$parent)
        kindata.means <- aggregate(kindata, list(kindata$t), mean)
	kindata.means.mean <- mean(kindata.means$parent, na.rm=TRUE)

	n.parms = length(coef(m))
	df = length(kindata.means$parent) - n.parms
	kindata.means$est <- predict(m, kindata.means)

	f <- function(err)
	{
		(sum((kindata.means$parent - kindata.means$est)^2/((err*kindata.means.mean)^2)) - 
		 qchisq(1 - alpha,df))^2
	}
	err.min <- optimize(f, c(0.01,0.9))$minimum
	return(err.min)
}
