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

"get.refrunno" <- function(database=".ref.db") {

  cat("Please enter the run number you want to compare the current\n")
  cat("model with (0 = exit): ")
  new.runno <- readline()
  if(new.runno == 0) return()
   
  ## Check if it exists
   
  if(exists(paste("xpdb", new.runno, sep = ""))) {
    cat("A matching database was found. Do you want to use it? (y) (r=recreate) ")
    ans <- readline()
    if(ans == "y" || ans == "") {
      ##
      ## Use old
      ##
      c1<-call("assign",pos = 1, database, eval(as.name(paste("xpdb", 
               new.runno, sep = ""))))
      eval(c1)
      return(new.runno)
    }
  }

  ##
  ## Recreate
  ##
  newdb <- xpose.data(new.runno)
  if(is.null(newdb)) {
    cat("No new database read.\n")
    return()
  } else {
    newnam <- paste("xpdb", new.runno, sep = "")
    c2<-call("assign",pos = 1, newnam, newdb)
    eval(c2)
    c3<-call("assign",pos = 1, database,newdb )
    eval(c3)
  }
  return(new.runno)
}
