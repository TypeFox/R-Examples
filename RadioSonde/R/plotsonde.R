plotsonde <- 
function (dataframe, skewT=TRUE, winds=FALSE, site = "", title = "", 
            windplot = NULL, s = 3., col = c(1, 2), ... ){
#
# Copyright 2001,2002 Tim Hoar, Eric Gilleland, and Doug Nychka
#
# This file is part of the RadioSonde library for R and related languages.
#
# RadioSonde is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# RadioSonde is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RadioSonde; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

    msg <- deparse(match.call())

    if( skewT & winds){

       #
       # Plot the SKEW-T, log p diagram and the wind profile.
       #

       # Need some room for both the skewT plot and the wind profile.

       mar.skewt <- c(5.0999999999999996, 1.1000000000000001, 
                      2.1000000000000001, 5.0999999999999996)
       skewt.plt <- skewt.axis(mar = mar.skewt)$plt
       title(title)

       if(is.null(windplot)) {
                windplot <- skewt.plt
                windplot[1] <- 0.8
                windplot[2] <- 0.95

       } else if( (windplot[2] < windplot[1]) | 
                  (windplot[4] < windplot[3]) ) {
                stop("plot region (windplot) too small to add second plot")
       }

       first.par <- par()

       # Draw the SKEW-T, log p diagram
       # Draw background and overlay profiles

       skewt.axis()
       skewt.lines(dataframe$temp,  dataframe$press, col = col[1], ...)
       skewt.lines(dataframe$dewpt, dataframe$press, col = col[2], ...)

       #
       # Draw the windplot in the "space allocated"
       # top and bottom mar the same as skewt
       #
        print( windplot)
        par(new = TRUE, pty = "m", plt = windplot, err = -1.)
        plotwind(dataframe = dataframe, size = s, legend = FALSE)
        par(plt = first.par$plt, mar = first.par$mar, new = FALSE, pty = first.par$
                pty, usr = first.par$usr)
        #       title1 <- paste(site, ": ", month.year, " ", time, sep = "")
        invisible()

    } else if( skewT & !winds) {

       #
       # Draw the SKEW-T, log p diagram
       # Draw background and overlay profiles
       #

       skewt.axis()
       skewt.lines(dataframe$temp,  dataframe$press, col = col[1], ...)
       skewt.lines(dataframe$dewpt, dataframe$press, col = col[2], ...)
       title(title)

    } else if( !skewT & winds) {

       #
       # Draw the Wind profile only
       #
       plotwind(dataframe=dataframe, ...)
       title(title)

    }  # end of if else stmts

    invisible()
}
