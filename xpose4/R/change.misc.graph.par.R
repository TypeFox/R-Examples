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

"change.misc.graph.par"  <- function(object, classic = FALSE)
{
  data <- object
  
  cat("These are the current general plot settings:\n\n")
  cat(paste("Type:",data@Prefs@Graph.prefs$type,sep=" "),"\n")
  cat(paste("Plotting character:", data@Prefs@Graph.prefs$pch,sep=" "),"\n")
  cat(paste("Point size:",data@Prefs@Graph.prefs$cex,sep=" "),"\n")
  cat(paste("Color:",data@Prefs@Graph.prefs$col,sep=" "),"\n")
  cat(paste("Line type:", data@Prefs@Graph.prefs$lty,sep=" "),"\n")
  cat(paste("Line width:",data@Prefs@Graph.prefs$lwd,sep=" "),"\n")
  cat(paste("Grid:",data@Prefs@Graph.prefs$grid,sep=" "),"\n")
  cat(paste("Aspect:",data@Prefs@Graph.prefs$aspect,sep=" "),"\n")
  cat("\n")

  # gr.stngs <- xp.gr.stngs

  cat("Specify a new plot type or leave blank to keep unchanged.\n\n")
  cat("(1-character string giving the type of plot desired.  The\n")
  cat("following values are possible, for details, see \'plot\': \'\"p\"\'\n")
  cat("for points, \'\"l\"\' for lines, \'\"b\"\' for both.): \n\n")
  
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$type <- ans
  }

  cat("\nSpecify a new plotting character or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer. See R help on \'points\'. The default is\n")
  cat("1, an open circle.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$pch <- as.numeric(ans)
  }

  cat("\nSpecify a new point scale or leave blank to keep unchanged.\n\n")
  cat("(The amount by which plotting text and symbols should be scaled\n")
  cat("relative to the default. 'NULL' and 'NA' are equivalent to '1.0'.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$cex <- as.numeric(ans)
  }
  
  cat("\nSpecify a new line color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is blue (col=4).): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$col <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$col <- ans
  }

  cat("\nSpecify a new line type or leave blank to keep unchanged.\n\n")
  cat("(Line types are specified as an integer (0=blank, 1=solid, \n")
  cat("2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).) \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$lty <- as.numeric(ans)
  }

  cat("\nSpecify a new line width or leave blank to keep unchanged.\n\n")
  cat("(A positive real number): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$lwd <- as.numeric(ans)
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
