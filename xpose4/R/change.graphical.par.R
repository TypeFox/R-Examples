# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"change.graphical.par"  <-
  function() {


    choices <- c("Return to the previous menu ->",
                 "Change general plot settings",
                 "Change settings for the smooth",
                 "Change settings for the linear regression line",
                 "Change settings for the line of identity",
                 "Change settings for labelling points",
                 "Change settings for box-and-whisker plots",
                 "Change dilution settings",
                 "Change prediction interval (PI) settings",
                 "Change conditioning settings",
                 "Reset graphical parameters to factory settings",
                 "Export graphical settings to a file",
                 "Import graphical settings from a file"
                 )


    title="\nGRAPHICAL PARAMETERS MENU\n  \\main\\preferences\\Change graphical parameters"
    
    pick <- menu(choices,title=title)
   
    qx <- 0
    switch(pick+1,
           qx <- 2,
           qx <- 1,
           change.misc.graph.par(get(".cur.db"), classic = T),
           change.smooth.graph.par(get(".cur.db"), classic = T),
           change.lm.graph.par(get(".cur.db"), classic = T),
           change.ab.graph.par(get(".cur.db"), classic = T),
           change.label.par(get(".cur.db"), classic = T),
           change.bw.graph.par(get(".cur.db"), classic = T),
           change.dil.graph.par(get(".cur.db"), classic = T),
           change.pi.graph.par(get(".cur.db"), classic = T),
           change.cond.graph.par(get(".cur.db"), classic = T),
           reset.graph.par(get(".cur.db"), classic = T),
           export.graph.par(get(".cur.db")),
           import.graph.par(get(".cur.db"), classic = T),
           )

    if(qx == 2) {
      return(invisible(2))
    } else {
      if(qx == 1) {
        return(invisible(0))
      } else {
        Recall()
      }
    } 

}
