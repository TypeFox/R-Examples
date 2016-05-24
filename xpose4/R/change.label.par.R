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

"change.label.par"  <- function(object, classic = FALSE)
{
  data <- object
  
  cat("These are the current labelling settings:\n\n")
  cat(paste("Label points:",data@Prefs@Graph.prefs$ids,sep=" "),"\n")
  cat(paste("Labelling style:", data@Prefs@Graph.prefs$idsmode,sep=" "),"\n")
  cat(paste("Extremes:",data@Prefs@Graph.prefs$idsext,sep=" "),"\n")
  cat(paste("Label size:",data@Prefs@Graph.prefs$idscex,sep=" "),"\n")
  cat(paste("Direction:", data@Prefs@Graph.prefs$idsdir,sep=" "),"\n")
  cat("\n")

  # gr.stngs <- xp.gr.stngs

  cat("Label the data points with ID values?\n\n")
  cat("(TRUE or NULL): \n\n")
  
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$ids <- ans
  }

  cat("Label all points, or extremes only?\n\n")
  cat("(NULL specifies extremes, and is the default. Anything else \n")
  cat("indicates that all points should be labelled.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$idsmode <- ans
  }

  cat("Specify the extremes to use for labelling points.\n\n")
  cat("(A number between 0 and 1. The default is 0.05 (only the most \n")
  cat("extreme 5% of points are labelled).): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$idsext <- as.numeric(ans)
  }
  
  cat("Specify a new point scale or leave blank to keep unchanged.\n\n")
  cat("(The amount by which plotting text and symbols should be scaled\n")
  cat("relative to the default. 'NULL' and 'NA' are equivalent to '1.0'.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$idscex <- as.numeric(ans)
  }

  cat("Specify the direction for labelling or leave blank to keep unchanged.\n\n")
  cat("(Possible values are 'up', 'down' and 'both'.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$idsdir <- ans
  }

  if (classic==TRUE) {
    c1<-call("assign",paste("xpdb", object@Runno, sep = ""), data, immediate=T, envir = .GlobalEnv)
    eval(c1)
    c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
    eval(c2)
    return(cat(""))
    
  } else {
    return(data)
  }

}
