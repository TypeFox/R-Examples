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
"tim.colors" <- function(n = 64, alpha = 1) {
    # tims original 64 color definition definition:
    orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
        "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
        "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
        "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
        "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
        "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
        "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
        "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
        "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
        "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
        "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
        "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
        "#AF0000", "#9F0000", "#8F0000", "#800000")
    if (n == 64 & alpha == 1) 
        return(orig)
    rgb.tim <- t(col2rgb(orig))
    temp <- matrix(NA, ncol = 3, nrow = n)
    x <- seq(0, 1, , 64)
    xg <- seq(0, 1, , n)
    for (k in 1:3) {
        hold <- splint(x, rgb.tim[, k], xg)
        hold[hold < 0] <- 0
        hold[hold > 255] <- 255
        temp[, k] <- round(hold)
    }
    if (alpha == 1) {
        rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    }
    else {
        rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255, 
            alpha = alpha)
    }
}
