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

"change.xlabel"<-
function(object, listall=TRUE, classic=FALSE)
{
	if(listall) {
		db.names(object)
	}
	
	cat("\nPlease type the name of the variable whose label you\n")
  cat("would like to change...\n")
  ans <- readline()
  dat <- object
        
  if(ans == "NULL" || ans == "null") {
	  dat@Prefs@Xvardef$idlab <- NULL
	  c1<-call("assign","object", dat, immediate=T, envir = .GlobalEnv)
    eval(c1)
	  c2<-call("assign",paste("xpdb", object@Runno, sep = ""), dat, immediate=T, envir = .GlobalEnv)
    eval(c2)
	  invisible()
	  return()
	  
	} else {
    
    if(is.na(pmatch(ans, names(dat@Data)))) {
		  cat("Couldn't find", ans, "in the current database!\n")
		  Recall(object, listall = F)
	  } else {
      if(ans == 0) {
		    return()
	    } else {
		    cat(paste("The current label for variable ",ans," is '",xlabel(ans, dat),"'.\n", sep=""))
		    cat("Enter the new one:\n")
		    labl <- readline()
		    dat@Prefs@Labels[ans] <- labl
		    
        if (classic==TRUE) {
          c1<-call("assign",paste("xpdb", object@Runno, sep = ""), data, immediate=T, envir = .GlobalEnv)
          eval(c1)
          c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
          eval(c2)
          return(cat(""))
        } else {
          return(dat)
        }
      }
	  }
  }
}
