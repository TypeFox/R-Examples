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

"add.dichot"  <- function(object, listall=TRUE, classic=FALSE )
{
  if(listall) db.names(object)
  
  data <- object@Data
  sdata <- object@SData
  
  nams <- names(data)
  
  cat("Please type the name of the item you whish to categorize:\n")
  ans <- readline()
  
  if(is.na(match(ans,nams))) {
    cat("Couldn't find that data item in the current database.\n")
    return(cat(""))
  }
  
  if(!is.numeric(data[[ans]])) {
    cat("The specified data item is not numeric.\n")
    return(cat(""))
  }

  cat("The quartiles of",ans,"are:\n")
  cat(paste(summary(data[[ans]])[1:3]),"\n\n",sep=" ")
  
  cat("Enter the breakpoints to use and finish with a blank line\n")
  cat("(the quartiles will be used if left empty).\n")
 
  br <- scan(what=numeric())

  if(length(br) == 0) {
    cat("No breakpoints were given. Will use default of Q1-3.\n")
    br <- summary(data[[ans]])[1:3]
  }
     
  br <- c(min(data[[ans]]),br,max(data[[ans]]))
  br <- unique(br)
  br <- br[order(br)]
  
  nam <- paste("cat", ans, sep="")
  data[[nam]] <- cut(data[[ans]], breaks=br, include.lowest=T)
  if (!is.null(sdata)) {
    sdata[[nam]] <- cut(sdata[[ans]], breaks=br, include.lowest=T)
  }
  
  object@Data <- data
  object@SData <- sdata
  
  #for (i in items) {
  catitem <- paste("cat", ans, sep="")
  object@Prefs@Labels[[catitem]] <- c(paste("categorical(", ans, ")", sep=""))
  #}

  if (classic==TRUE) {
    #     assign(paste("xpdb", object@Runno, sep = ""), object, immediate=T, envir = .GlobalEnv)
    #     assign(pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
    #     return(cat(""))
    
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

