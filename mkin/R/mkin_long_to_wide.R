# $Id$

# Copyright (C) 2010-2011 Johannes Ranke
# Contact: mkin-devel@lists.berlios.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

mkin_long_to_wide <- function(long_data, time = "time", outtime = "time")
{
  colnames <- unique(long_data$name)
  wide_data <- data.frame(time = subset(long_data, name == colnames[1], time))
  names(wide_data) <- outtime
  for (var in colnames) {
    wide_data[var] <- subset(long_data, name == var, value)
  }
  return(wide_data)
}
