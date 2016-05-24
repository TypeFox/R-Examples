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

"change.pi.graph.par"  <- function(object, classic = FALSE)
{
  data <- object
  
  cat("These are the current settings for prediction intervals:\n\n")
  cat(paste("PI limits:",data@Prefs@Graph.prefs$PIlimits,sep=" "),"\n")
  cat(paste("Upper limit plot type:",data@Prefs@Graph.prefs$PIuptyp,sep=" "),"\n")
  cat(paste("Lower limit plot type:",data@Prefs@Graph.prefs$PIdotyp,sep=" "),"\n")
  cat(paste("Median plot type:",data@Prefs@Graph.prefs$PImetyp,sep=" "),"\n")
  cat(paste("Upper limit line colour:",data@Prefs@Graph.prefs$PIupcol,sep=" "),"\n")
  cat(paste("Lower limit line colour:",data@Prefs@Graph.prefs$PIdocol,sep=" "),"\n")
  cat(paste("Median line colour:",data@Prefs@Graph.prefs$PImecol,sep=" "),"\n")
  cat(paste("Upper limit line type:",data@Prefs@Graph.prefs$PIuplty,sep=" "),"\n")
  cat(paste("Lower limit line type:",data@Prefs@Graph.prefs$PIdolty,sep=" "),"\n")
  cat(paste("Median line type:",data@Prefs@Graph.prefs$PImelty,sep=" "),"\n")
  cat(paste("Upper limit line width:",data@Prefs@Graph.prefs$PIuplwd,sep=" "),"\n")
  cat(paste("Lower limit line width:",data@Prefs@Graph.prefs$PIdolwd,sep=" "),"\n")
  cat(paste("Median line width:",data@Prefs@Graph.prefs$PImelwd,sep=" "),"\n")
  cat("\n")

  # gr.stngs <- xp.gr.stngs
  
  cat("Specify the upper limit of the prediction interval.\n\n")
  cat("(a number, between 0 and 1. The default is 0.975.): \n\n")
  ans <- readline()
  if(ans!="") {
    PIuplimit <- as.numeric(ans)
  }
  
  cat("\nSpecify the lower limit of the prediction interval.\n\n")
  cat("(a number, between 0 and 1. The default is 0.025.): \n\n")
  ans <- readline()
  if(ans!="") {
    PIdolimit <- as.numeric(ans)
  }
  
  data@Prefs@Graph.prefs$PIlimits = c(PIdolimit, PIuplimit) 
  
  cat("\nSpecify a new plot type for the upper PI limit, or leave blank to keep\n")
  cat("unchanged.\n\n")
  cat("(1-character string giving the type of plot desired.  The\n")
  cat("following values are possible, for details, see \'plot\': \'\"p\"\'\n")
  cat("for points, \'\"l\"\' for lines, or \'\"b\"\' for both.): \n\n")
  
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PIuptyp <- ans
  }

  cat("\nSpecify a new plot type for the lower PI limit, or leave blank to keep\n")
  cat("unchanged.\n\n")
  cat("(1-character string giving the type of plot desired.  The\n")
  cat("following values are possible, for details, see \'plot\': \'\"p\"\'\n")
  cat("for points, \'\"l\"\' for lines, or \'\"b\"\' for both.): \n\n")
  
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PIdotyp <- ans
  }
  
  cat("\nSpecify a new plot type for the median, or leave blank to keep\n")
  cat("unchanged.\n\n")
  cat("(1-character string giving the type of plot desired.  The\n")
  cat("following values are possible, for details, see \'plot\': \'\"p\"\'\n")
  cat("for points, \'\"l\"\' for lines, or \'\"b\"\' for both.): \n\n")
  
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PImetyp <- ans
  }
  
  cat("\nSpecify a new upper limit line color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is black.): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$PIupcol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$PIupcol <- ans
  }
  
  cat("\nSpecify a new lower limit line color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is black.): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$PIdocol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$PIdocol <- ans
  }
  
  cat("\nSpecify a new median line color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is black.): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$PImecol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$PImecol <- ans
  }
  
  cat("\nSpecify a new upper limit line type or leave blank to keep unchanged.\n\n")
  cat("(Line types are specified as an integer (0=blank, 1=solid, \n")
  cat("2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).) \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PIuplty <- as.numeric(ans)
  }
  
  cat("\nSpecify a new lower limit line type or leave blank to keep unchanged.\n\n")
  cat("(Line types are specified as an integer (0=blank, 1=solid, \n")
  cat("2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).) \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PIdolty <- as.numeric(ans)
  }
  
  cat("\nSpecify a new median line type or leave blank to keep unchanged.\n\n")
  cat("(Line types are specified as an integer (0=blank, 1=solid, \n")
  cat("2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).) \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PImelty <- as.numeric(ans)
  }

  cat("\nSpecify a new upper limit line width or leave blank to keep unchanged.\n\n")
  cat("(A positive real number): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PIuplwd <- as.numeric(ans)
  }
  
  cat("\nSpecify a new lower limit line width or leave blank to keep unchanged.\n\n")
  cat("(A positive real number): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PIdolwd <- as.numeric(ans)
  }
  
  cat("\nSpecify a new median line width or leave blank to keep unchanged.\n\n")
  cat("(A positive real number): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$PImelwd <- as.numeric(ans)
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
