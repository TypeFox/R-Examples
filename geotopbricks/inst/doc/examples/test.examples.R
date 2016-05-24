# file test.examples.R
#
# This file contains a script which plots the measurement sites on a GoogleMap support
#
#
# author: Emanuele Cordano on 08-11-2012

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################



library(geotopbricks)


dir <- system.file("doc/examples",package="geotopbricks")
this_file <- system.file("doc/examples/test.examples.R",package="geotopbricks")

l <- list.files(dir,full.name=TRUE)
l <- l[l!=this_file]

for (it in l) source(it)
