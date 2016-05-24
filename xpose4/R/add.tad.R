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

"add.tad"<-
  #function(runno = get(fr = 0, "runno"))
  function(object, classic=FALSE)
{

  tab.suffix <- ""
  ans <- ""

  idv.nam <- object@Prefs@Xvardef$idv
  if(is.null(idv.nam)) {
    cat("The IDV wasn't set in the current database.\n")
    cat("Please type the name of the independent variable to use in\n")
    cat("calculations, or press return for the default table file column\n")
    cat("(sdtab and mytab column 2, mutab column 3): ")

    ans <- readline()
  }
  
  "make.tad" <- function(indat, idvcol) {
    tad <- rep(0, length.out = nrow(indat))
    for(i in 1:nrow(indat)) {
      if(indat$WRES[i] == 0) {
        last.dose <- indat[[idvcol]][i]
      } else {
        tad[i] <- indat[[idvcol]][i] - last.dose
      }
    }
    return(tad)
  }

  dat <- object@Data
  sdat <- object@SData

  dat$TAD <- make.tad(dat, match(xvardef("idv", object), names(dat)))
  if(!is.null(sdat))
    sdat$TAD <- make.tad(sdat, match(xvardef("idv", object), names(sdat)))
  
  object@Data <- dat

  if(!is.null(sdat))
    object@SData <- sdat
  
  data <- object
  
  data@Prefs@Labels$TAD <- c("Time after dose (h)")

  #assign(object, data, immediate=T, env = .GlobalEnv)
  
  if (classic==TRUE) {
    #assign(paste("xpdb", object@Runno, sep = ""), data, immediate=T, envir = .GlobalEnv)
    #assign(pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
    #return(cat(""))
    
    ## to avoid checks on global variable assignment in package building
    c1<-call("assign",paste("xpdb", object@Runno, sep = ""),object,envir=.GlobalEnv)
    eval(c1)
    c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
    eval(c2)
    return(cat(""))
    
  } else {
    return(data)
  }
}
