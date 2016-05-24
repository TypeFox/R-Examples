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

"change.ab.graph.par"  <- function(object, classic = FALSE)
{
  data <- object
  
  cat("These are the current settings for lines of identity:\n\n")
  cat(paste("Use line of identity:",data@Prefs@Graph.prefs$abline,sep=" "),"\n")
  cat(paste("Line color:",data@Prefs@Graph.prefs$ablcol,sep=" "),"\n")
  cat(paste("Line type:", data@Prefs@Graph.prefs$abllty,sep=" "),"\n")
  cat(paste("Line width:",data@Prefs@Graph.prefs$abllwd,sep=" "),"\n")
  cat("\n")

  # gr.stngs <- xp.gr.stngs
  
  cat("Should a line of identity be included by default?\n\n")
  cat("(A 2-element list specifying whether, and if so, what\n")
  cat("type of line to automatically overlay on plots. NULL indicates\n")
  cat("no line unless overridden by a specific function (e.g. dv.vs.pred),\n")
  cat("c(1,1) indicates an x:y line, and c(1,2) indicates an x:2y line,\n")
  cat("for example. See R help on 'panel.abline' for more information.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$abline <- ans
  }

  cat("\nSpecify a new line color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is black (col=1).): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$ablcol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$ablcol <- ans
  }

  cat("\nSpecify a new line type or leave blank to keep unchanged.\n\n")
  cat("(Line types are specified as an integer (0=blank, 1=solid, \n")
  cat("2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).) \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$abllty <- as.numeric(ans)
  }

  cat("\nSpecify a new line width or leave blank to keep unchanged.\n\n")
  cat("(A positive real number): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$abllwd <- as.numeric(ans)
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
