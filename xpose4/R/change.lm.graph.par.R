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

"change.lm.graph.par"  <- function(object, classic = FALSE)
{
  data <- object
  
  cat("These are the current settings for linear regression lines:\n\n")
  cat(paste("Use linear regression line:",data@Prefs@Graph.prefs$lmline,sep=" "),"\n")
  cat(paste("Line color:",data@Prefs@Graph.prefs$lmcol,sep=" "),"\n")
  cat(paste("Line type:", data@Prefs@Graph.prefs$lmlty,sep=" "),"\n")
  cat(paste("Line width:",data@Prefs@Graph.prefs$lmlwd,sep=" "),"\n")
  cat("\n")

  # gr.stngs <- xp.gr.stngs
  
  cat("Should linear regression lines be included by default?\n\n")
  cat("(TRUE or FALSE(=NULL).): \n\n")
  ans <- readline()
  if(ans!="") {
    if (ans=="NULL") {
      ans = "FALSE"
    }
    data@Prefs@Graph.prefs$suline <- as.logical(ans)
  }

  cat("\nSpecify a new line color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is blue.): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$lmcol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$lmcol <- ans
  }

  cat("\nSpecify a new line type or leave blank to keep unchanged.\n\n")
  cat("(Line types are specified as an integer (0=blank, 1=solid, \n")
  cat("2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).) \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$lmlty <- as.numeric(ans)
  }

  cat("\nSpecify a new line width or leave blank to keep unchanged.\n\n")
  cat("(A positive real number): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$lmlwd <- as.numeric(ans)
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
