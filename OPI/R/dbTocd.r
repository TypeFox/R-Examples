# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: June 2012
#
# Modified Tue  8 Jul 2014: make generic for any dB scale
#
# Copyright 2012 Andrew Turpin and Jonathan Denniss
# This program is part of the OPI (http://perimetry.org/OPI).
# OPI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

    # Convert dB to cd/m^2
    # Default is HFA units (maxStim = 10000 apostilbs)
    # maxStim is in cd/m^2
dbTocd <- function(db, maxStim=10000/pi) { maxStim * 10^(-db/10) }

    # Convert cd/m^2 to dB
    # default is HFA units (maxStim = 10000 apostilbs)
    # maxStim is in cd/m^2
cdTodb <- function(cd, maxStim=10000/pi) { -10*log10(cd/maxStim) }
