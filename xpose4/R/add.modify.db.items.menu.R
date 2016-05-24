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

"add.modify.db.items.menu" <-
  function() {
    
    choices <- c("Return to the previous menu ->",
                 "Add a time after dose (TAD) item to the current database", 
                 "*Add averaged covariates to the current database", 
                 "Add exponentiated values of an item to the current database",
                 "Add the logarithm of an item to the current database",
                 "Add absolute values of an item to the current database",
                 "* Add a random bivariate factor to the current database",
                 "Change covariates categorical <-> continuous ", 
                 "Create a categorical item based on a continuous item"
                 ##"Copy data item from another database"
    )
    
    title=("\nCREATE/MODIFY DATA ITEMS MENU\n  \\main\\preferences\\Create/modify items in the current database")
    pick <- menu(choices,title=title)
    qx <- 0
    switch(pick+1,
           qx <- 2,
           qx <- 1,
           add.tad(eval(parse(text=".cur.db")), classic=T),
           cat("Not yet implemented!\n"), #cov.ave(),
           add.exp(eval(parse(text=".cur.db")), classic=T),
           add.log(eval(parse(text=".cur.db")), classic=T),
           add.absval(eval(parse(text=".cur.db")), classic=T),
           #object<-add.absval(object),
           cat("Not yet implemented!\n"), #add.rand(get("eval(parse(text=".cur.db"))", classic=T),
           change.cat.cont(eval(parse(text=".cur.db")), classic=T),
           add.dichot(eval(parse(text=".cur.db")), classic=T)
           ##copy.item()
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

