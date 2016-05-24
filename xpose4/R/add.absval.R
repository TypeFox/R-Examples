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

"add.absval" <- function(object, listall=TRUE, classic=FALSE )
{
  if(listall) db.names(object)

  cat("Please type the names of the items to be converted to\n")
  cat("absolute values, one per line, and finish with a blank line.\n")

  items <- scan(what=character())
  data <- object@Data
  sdata <- object@SData
  nams <- names(data)
  for(i in items) {
    if(is.na(match(i, nams))) {
      cat("No match: ", i, "\n", sep="")
      next
    }

    nam <- paste("abs", i, sep="")
    #cat(nam)
    newit <- abs(data[[i]])
    if(any((data[[i]]-newit) !=0)) {
      data[[nam]] <- newit
      #vname(data[,nam]) <- nam
    }
    if(!all(is.null(sdata))){ 
       newits <- abs(sdata[[i]])
       if(any((sdata[[i]]-newits) !=0)) {
         sdata[[nam]] <- newits
         #vname(data[,nam]) <- nam
       }
    }
  }
  object@Data <- data
  object@SData <- sdata
  
  for (i in items) {
    expitem <- paste("abs", i, sep="")
    object@Prefs@Labels[[expitem]] <- c(paste("abs(", i, ")", sep=""))
  }
  
  if (classic==TRUE) {
    c1<-call("assign",paste("xpdb", object@Runno, sep = ""),object,envir=.GlobalEnv)
    eval(c1)
    #assign(paste("xpdb", object@Runno, sep = ""), object, immediate=T, envir = .GlobalEnv)
    c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
    eval(c2)
    #assign(pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))    
    return(cat(""))
  } else {
    return(object)
  }

}
