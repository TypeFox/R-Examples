#  Part of the devEMF Package for R.  Copyright 2011 by Philip Johnson.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

emf <- function(file = "Rplot.emf", width=7, height=7,
                bg = "transparent", fg = "black", pointsize=12,
                family = "Helvetica", custom.lty=FALSE)
{
  ps.options() #initialize fonts for retrieving metric information
  .External('devEMF', PACKAGE='devEMF', file, bg, fg, width, height, pointsize,
            family, custom.lty)
  invisible()
}
