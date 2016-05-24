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

"main.menu" <- function() {

  choices <- c("Documentation ->",
               "Preferences ->",
               "Data checkout ->",
               "Goodness of fit plots ->",
               "Parameters ->",
               "Covariate model ->",
               "Model comparison ->",
               "Conditional weighted residuals ->",
               "Visual and numerical predictive check plots ->",
               "License and citation information",
               "Quit")

  title=paste(
    "\nMAIN MENU",
    "\n  Enter an item from the menu, or 0 to exit",
    "\n  -> : Indicates a directory",
    "\n  *  : Indicates functionality not yet available",
    sep="")

  pick <- menu(choices,title=title)

  qx <- 0
  switch(pick + 1,
         qx <- 1,
         qx <- documentation.menu(),
         qx <- preferences.menu(),
         qx <- data.checkout.menu(),
         qx <- gof.menu(),
         qx <- parameters.menu(),
         qx <- covariate.model.menu(),
         qx <- model.comparison.menu(),
         qx <- cwres.menu(),
         qx <- vpc.npc.menu(),
         xpose.license.citation(),
         qx <- 1
         )
  

  # quit menu system
  if(qx == 1 || qx==2) {
    ## Turn of graphics window
    ##if(dev.cur() > 1){
    ##  dev.off()
    ##}
    return(invisible(1))
  } else {
    Recall()
  }
  
}
