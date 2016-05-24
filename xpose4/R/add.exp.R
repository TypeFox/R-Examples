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

"add.exp" <- function(object, listall=TRUE, classic=FALSE )
{
  if(listall) db.names(object)

  cat("Please type the names of the items to be exponentiated, one\n")
  cat("per line, and finish with a blank line.\n")

  items <- scan(what=character())
  data <- object@Data
  sdata <- object@SData
  nams <- names(data)
  for(i in items) {
    if(is.na(match(i, nams))) {
      cat("No match: ", i, "\n", sep="")
      next
    }

    nam <- paste("exp", i, sep="")
    #cat(nam)
    newit <- exp(data[[i]])
    newits <- exp(sdata[[i]])
    if(any((data[[i]]-newit) !=0)) {
      data[[nam]] <- newit
      #vname(data[,nam]) <- nam
    }
    if(any((sdata[[i]]-newits) !=0)) {
      sdata[[nam]] <- newits
      #vname(data[,nam]) <- nam
    }
  }

  object@Data <- data
  object@SData <- sdata
  
  for (i in items) {
    expitem <- paste("exp", i, sep="")
    object@Prefs@Labels[[expitem]] <- c(paste("exp(", i, ")", sep=""))
  }
  
  if (classic==TRUE) {
    #assign(paste("xpdb", object@Runno, sep = ""), object, immediate=T, envir = .GlobalEnv)
    #assign(pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
    #return(cat(""))
    
    ## to avoid checks on global variable assignment in package building
    c1<-call("assign",paste("xpdb", object@Runno, sep = ""),object,envir=.GlobalEnv)
    eval(c1)
    c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
    eval(c2)
    return(cat(""))
    
  } else {
    return(object)
  }

}
