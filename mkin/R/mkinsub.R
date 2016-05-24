# Copyright (C) 2014,2015 Johannes Ranke
# Portions of this code are copyright (C) 2013 Eurofins Regulatory AG
# Contact: jranke@uni-bremen.de

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
mkinsub <- function(submodel, to = NULL, sink = TRUE, full_name = NA)
{
  return(list(type = submodel, to = to, sink = sink, full_name = full_name))
}
# vim: set ts=2 sw=2 expandtab:
