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

"change.dil.graph.par"  <- function(object, classic = FALSE)
{
  data <- object
  
  cat("These are the current data dilution settings:\n\n")
  cat(paste("Dilution type:",data@Prefs@Graph.prefs$diltype,sep=" "),"\n")
  cat(paste("Fraction to be diluted:",data@Prefs@Graph.prefs$dilfrac,sep=" "),"\n")
  cat(paste("Dilution confidence interval:", data@Prefs@Graph.prefs$dilci,sep=" "),"\n")
  cat("\n")

  # gr.stngs <- xp.gr.stngs
  
  cat("Stratified dilution?\n\n")
  cat("(NULL specifies random dilution without stratification, and is the\n")
  cat("default. Anything else indicates that stratified dilution is to be used.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$diltype <- ans
  }
  
  cat("Specify the fraction eligible for dilution.\n\n")
  cat("(a number, between 0 and 1, defining the fraction of the data to\n") 
  cat("be diluted. The default is 0.7.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$dilfrac <- as.numeric(ans)
  }
  
  cat("Specify the confidence interval for stratified dilution.\n\n")
  cat("(a number, between 0 and 1, defining the range of the data\n") 
  cat("eligible for dilution. The default is 0.95.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$dilci <- as.numeric(ans)
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
