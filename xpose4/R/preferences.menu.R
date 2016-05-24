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

preferences.menu <-
  function() {

    choices <- c("Return to the main menu ->",
                 "Change run number/database", 
                 "Manage variables in the current database ->",
                 "Create/modify items in the current database ->",
                                        #"Variable preferences ->",
                 "Change graphical parameters (line colors, widths, etc) ->"
                 )

    title="\nPREFERENCES MENU\n  \\main\\preferences"
    
    pick <- menu(choices,title=title)
   
    qx <- 0
    switch(pick+1,
           qx <- 2,
           qx <- 1,
           change.xp.obj(),
           qx <- manage.db(),
           qx <- add.modify.db.items.menu(),
                                        #qx <- table.file.read.settings.menu(),
           qx <- change.graphical.par()
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
