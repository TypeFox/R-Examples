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

"change.cond.graph.par"  <- function(object, classic = FALSE)
{
  data <- object
  
  cat("These are the current conditioning settings:\n\n")
  cat(paste("Condition on:",data@Prefs@Graph.prefs$condvar,sep=" "),"\n")
  cat(paste("Order by:",data@Prefs@Graph.prefs$ordby,sep=" "),"\n")
  cat(paste("Ordering function:",data@Prefs@Graph.prefs$byordfun,sep=" "),"\n")
  cat(paste("Default number of shingles:", data@Prefs@Graph.prefs$shingnum,sep=" "),"\n")
  cat(paste("Shingle overlap:", data@Prefs@Graph.prefs$shingol,sep=" "),"\n")
  cat("\n")

  # gr.stngs <- xp.gr.stngs
  
  cat("Use a variable to condition plots made by Xpose?\n\n")
  cat("(A string of the variable name to condition on, can also\n")
  cat("be NULL, press ENTER/RETURN to leave it as it is): \n\n") 
  ans <- readline()
  if(ans == "NULL") {
    ans <- NULL
    data@Prefs@Graph.prefs$condvar <- ans
  } else {
    if(ans!="") {
      data@Prefs@Graph.prefs$condvar <- ans
    }
  }

  cat("Use a variable to reorder factorial variables during conditioning?\n\n")
  cat("(A string with the name of a variable to be used to reorder\n") 
  cat("any factorial conditioning variables. The variable is used in a\n") 
  cat("call to the R 'reorder.factor' function.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$ordby <- ans
  }
  
  cat("Specify a function for use with conditioning using categorical\n")
  cat("variables.\n\n")
  cat("(The name of the function to be used when reordering a\n") 
  cat("factor conditioning variable. Can be 'mean', 'median', etc.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$byordfun <- ans
  }
  
  cat("Specify the number of shingles to be used when conditioning on a\n")
  cat("continuous variable.\n\n")
  cat("(The number of shingles ('parts') a continuous\n")
  cat("conditioning variable should be divided into.): \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$shingnum <- as.numeric(ans)
  }
  
  cat("Specify the amount of overlap between adjacent shingles: \n\n")
  ans <- readline()
  if(ans!="") {
    data@Prefs@Graph.prefs$shingol <- as.numeric(ans)
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
