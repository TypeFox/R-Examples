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

"change.steppit"<-
function(first=TRUE)
{
  value <- .cur.db@Prefs@Gam.prefs$steppit
  if (first==TRUE){
    cat("Use a stepwise search (true/false)\n")
    if (is.null(value)) {
      cat("The current value is NULL...\n")
    } else {
      cat("The current value is",value,"...\n")
    }
    cat("\nPlease type the new value ")
  }
  
  ans <- readline()
        
  if(ans == "FALSE" || ans == "FALSE") {
    .cur.db@Prefs@Gam.prefs$steppit <- FALSE
    c1<-call("assign",pos = 1, ".cur.db", .cur.db)
    eval(c1)
    invisible()
    return()
    
  } else {
    if(ans == "true" || ans == "TRUE") {
      .cur.db@Prefs@Gam.prefs$steppit <- TRUE
      c1<-call("assign",pos = 1, ".cur.db", .cur.db)
      eval(c1)
      invisible()
      return()
    } else {
      cat("Please enter TRUE or FALSE ")
      Recall(first=FALSE)
    }
  }
}
