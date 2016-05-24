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

"change.bw.graph.par"  <- function(object, classic = FALSE)
{
  data <- object
  
  cat("These are the current settings for the box-and-whisker plots:\n\n")
  cat(paste("Horizontal:",data@Prefs@Graph.prefs$bwhoriz,sep=" "),"\n")
  cat(paste("Height-to-interbox space ratio:",data@Prefs@Graph.prefs$bwratio,sep=" "),"\n")
  cat(paste("Proportional widths:", data@Prefs@Graph.prefs$bwvarwid,sep=" "),"\n")
  cat(paste("Dot plotting character:",data@Prefs@Graph.prefs$bwdotpch,sep=" "),"\n")
  cat(paste("Dot colour:",data@Prefs@Graph.prefs$bwdotcol,sep=" "),"\n")
  cat(paste("Dot scale:",data@Prefs@Graph.prefs$bwdotcex,sep=" "),"\n")
  cat(paste("Rectangle fill colour:",data@Prefs@Graph.prefs$bwrecfill,sep=" "),"\n")
  cat(paste("Rectangle line colour:",data@Prefs@Graph.prefs$bwreccol,sep=" "),"\n")
  cat(paste("Rectangle line type:",data@Prefs@Graph.prefs$bwreclty,sep=" "),"\n")
  cat(paste("Rectangle line width:", data@Prefs@Graph.prefs$bwreclwd,sep=" "),"\n")
  cat(paste("Umbrella line colour:",data@Prefs@Graph.prefs$bwumbcol,sep=" "),"\n")
  cat(paste("Umbrella line type:",data@Prefs@Graph.prefs$bwumblty,sep=" "),"\n")
  cat(paste("Umbrella line width:",data@Prefs@Graph.prefs$bwumblwd,sep=" "),"\n")
  cat(paste("Outlier plotting character:",data@Prefs@Graph.prefs$bwoutpch,sep=" "),"\n")
  cat(paste("Outlier colour:", data@Prefs@Graph.prefs$bwoutcol,sep=" "),"\n")
  cat(paste("Outlier scale:",data@Prefs@Graph.prefs$bwoutcex,sep=" "),"\n")
  cat("\n")

  # gr.stngs <- xp.gr.stngs
  
  cat("Should the plots be horizontal?\n\n")
  cat("(TRUE or FALSE.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwhoriz <- as.logical(ans)
  }

  cat("\nSpecify a new box height to interbox space ratio, or leave blank\n")
  cat("to keep unchanged.\n\n")
  cat("(Specified as number. The default is 1.5.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwratio <- as.numeric(ans)
  }
  
  cat("Should the box widths be proportional?\n\n")
  cat("(If TRUE, widths of boxplots are proportional to the number of points\n")
  cat("used in creating them. The default is FALSE.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwvarwid <- as.logical(ans)
  }
  
  cat("\nSpecify a new plotting character for dots or leave blank to keep.\n\n")
  cat("unchanged.\n\n")
  cat("(Specified as an integer. See R help on \'points\'. The default is\n")
  cat("16.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwdotpch <- as.numeric(ans)
  }
  
  cat("\nSpecify a new dot color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is black (col=1).): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$bwdotcol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$bwdotcol <- ans
  }
  
  cat("\nSpecify a new dot scale or leave blank to keep unchanged.\n\n")
  cat("(The amount by which plotting text and symbols should be scaled\n")
  cat("relative to the default. 'NULL' and 'NA' are equivalent to '1.0'.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwdotcex <- as.numeric(ans)
  }

  cat("\nSpecify a new rectangle fill color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is transparent.): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$bwrecfill <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$bwrecfill <- ans
  }

  cat("\nSpecify a new rectangle line color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is blue (col=2).): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$bwreccol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$bwreccol <- ans
  }

  cat("\nSpecify a new rectangle line type or leave blank to keep unchanged.\n\n")
  cat("(Line types are specified as an integer (0=blank, 1=solid, \n")
  cat("2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).) \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwreclty <- as.numeric(ans)
  }

  cat("\nSpecify a new rectangle line width or leave blank to keep unchanged.\n\n")
  cat("(A positive real number): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwreclwd <- as.numeric(ans)
  }
  
  cat("\nSpecify a new umbrella line color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is blue (col=2).): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$bwumbcol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$bwumbcol <- ans
  }

  cat("\nSpecify a new umbrella line type or leave blank to keep unchanged.\n\n")
  cat("(Line types are specified as an integer (0=blank, 1=solid, \n")
  cat("2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash).) \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwumblty <- as.numeric(ans)
  }

  cat("\nSpecify a new umbrella line width or leave blank to keep unchanged.\n\n")
  cat("(A positive real number): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwumblwd <- as.numeric(ans)
  }
  
  cat("\nSpecify a new plotting character for outliers or leave blank to keep.\n\n")
  cat("unchanged.\n\n")
  cat("(Specified as an integer. See R help on \'points\'. The default is\n")
  cat("1.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwdotpch <- as.numeric(ans)
  }
  
  cat("\nSpecify a new outlier color or leave blank to keep unchanged.\n\n")
  cat("(Specified as an integer or a text string. A full list is obtained \n")
  cat("by the R command 'colours()'. The default is blue.): \n\n")
  ans <- readline()
  if ((ans!="") && (!is.na(as.numeric(ans)))) {
    data@Prefs@Graph.prefs$bwdotcol <- as.numeric(ans)
  } else {
    data@Prefs@Graph.prefs$bwdotcol <- ans
  }
  
  cat("\nSpecify a new outlier scale or leave blank to keep unchanged.\n\n")
  cat("(The amount by which plotting text and symbols should be scaled\n")
  cat("relative to the default. 'NULL' and 'NA' are equivalent to '1.0'.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$bwdotcex <- as.numeric(ans)
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
