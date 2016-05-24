# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"world" <- function(...) {
    map("world", ...)
    invisible()
}
world.color <- function(...) {
    cat("world.color has been depreciated. Please use fill options in\nthe world/map function.", 
        fill = TRUE)
}
world.land <- function(...) {
    cat("world.land has been depreciated. Please use fill options in\nthe world/map function.", 
        fill = TRUE)
}

in.land.grid <- function(...) {
    cat("world.land has been depreciated. Please refer to fields 6.7.1 or earlier to acces this function.", 
        fill = TRUE)
}


