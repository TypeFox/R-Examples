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

"gof.menu" <- function() {

  choices <- c("Return to previous menu ->",
               "Basic goodness of fit plots",
               "Dependent variables vs predictions",
               "Predictions vs independent variable",
               "Individual plots",
               "Structural model diagnostics ->",
               "Residual error model diagnostics ->"
               )

  title="\nGOODNESS OF FIT PLOTS MENU\n  \\main\\goodness of fit plots"

  pick <- menu(choices,title=title)

  qx <- 0
  switch(pick+1,
         qx <- 2,
         qx <- 1,
         print(basic.gof(eval(parse(text=".cur.db")))),
         print(dv.vs.pred.ipred(eval(parse(text=".cur.db")))),
         print(dv.preds.vs.idv(eval(parse(text=".cur.db")))),
         print(ind.plots(eval(parse(text=".cur.db")))),
         qx <- structural.diagnostics.menu(),
         qx <- residual.diagnostics.menu()
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
