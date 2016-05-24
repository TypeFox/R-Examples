# $Id: test_FOCUS_2006_A.R 59 2010-07-28 12:29:15Z jranke $

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

library(kinfit)
data(FOCUS_2006_A)
fits <- kinfit(FOCUS_2006_A, kinmodels = c("SFO", "HS"))
print(kinresults(fits)$results, digits=5)
print(kinresults(fits)$stats, digits=5)

data(FOCUS_2006_B)
fits <- kinfit(FOCUS_2006_B, kinmodels = c("SFO", "FOMC", "DFOP"))
print(kinresults(fits)$results, digits=5)
print(kinresults(fits)$stats, digits=5)

data(FOCUS_2006_C)
fits <- kinfit(FOCUS_2006_C, kinmodels = c("SFO", "FOMC", "DFOP"))
print(kinresults(fits)$results, digits=5)
print(kinresults(fits)$stats, digits=5)

data(FOCUS_2006_D)
fits <- kinfit(FOCUS_2006_D, kinmodels = c("SFO", "FOMC"))
print(kinresults(fits)$results, digits=5)
print(kinresults(fits)$stats, digits=3)

data(FOCUS_2006_E)
fits <- kinfit(FOCUS_2006_E, kinmodels = c("SFO", "FOMC", "DFOP"))
print(kinresults(fits)$results, digits=5)
print(kinresults(fits)$stats, digits=5)
