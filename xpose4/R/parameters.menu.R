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

parameters.menu <-
  function() {

    choices <- c("Return to previous menu ->",
                 "Numerically summarize the parameters",
                 "Distribution of parameters (QQ plots)",
                 "Distribution of parameters (histograms)",
                 "Scatterplot matrix of parameters",
                 "Parameter vs parameter",
                 "Distribution of random parameters (QQ plots)",
                 "Distribution of random parameters (histograms)",
                 "Scatterplot matrix of random parameters",
                 "* Random effects vs typical parameter values"
                 )

    title="\nPARAMETER PLOTS MENU\n  \\main\\parameters"

    pick <- menu(choices,title=title)
    
    qx <- 0
    switch(pick+1,
           qx <- 2,
           qx <- 1,
           parm.summary(eval(parse(text=".cur.db")),out.file=".ask"),
           print(parm.qq(eval(parse(text=".cur.db")),prompt=FALSE)),
           print(parm.hist(eval(parse(text=".cur.db")),prompt=FALSE)),
           print(parm.splom(eval(parse(text=".cur.db")))),
           print(parm.vs.parm(eval(parse(text=".cur.db")))),
           print(ranpar.qq(eval(parse(text=".cur.db")))),
           print(ranpar.hist(eval(parse(text=".cur.db")))),
           print(ranpar.splom(eval(parse(text=".cur.db")))),
           cat("\nNot defined yet\n")
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
