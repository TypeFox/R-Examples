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

runsum.print <- function(object,
                         ##dev="win",
                         ##printer="HP8150",
                         modfile=paste("run",object@Runno,".mod",sep=""),
                         listfile=paste("run",object@Runno,".lst",sep=""),
                         ##new.version=TRUE,
                         print.cex=0.45,
                         print.columns=3,
                         ...) {


  ## Get the name of the model file to use
  cat("Type the name of the model file (0=cancel, return=",modfile,")\n",sep="")
  ans <- readline()

  cmdfile <- NULL
  if(ans==0) {
    return()
  } else if (ans=="") {
    if(is.readable.file(modfile)) {
      cmdfile <- modfile
    }
  } else {
    if(is.readable.file(ans)) {
      cmdfile <- ans
    }
  }
  
  if(is.null(cmdfile)) {
    cat("The specified file couldn't be found in the current directory.\n")
    return()
  }

  ## Get the name of the list file to use
  cat("Type the name of the output file (0=cancel, return=",listfile,")\n",sep="")
  ans <- readline()

  lstfile <- NULL
  if(ans==0) {
    return()
  } else if (ans=="") {
    if(is.readable.file(listfile)) {
      lstfile <- listfile
    }
  } else {
    if(is.readable.file(ans)) {
      lstfile <- ans
    }
  }
  
  if(is.null(lstfile)) {
    cat("The specified file couldn't be found in the current directory.\n")
    return()
  }


  ## Start the printing device
  cat("Do you want to optimize the summary output for printing n(y)?\n")
  printit <- readline()
  if(printit=="y") { 
      dev.new(height=11.7,width=8.25)
  }

  if(printit=="y") { 
      runsum(object,
             modfile=cmdfile,
             listfile = lstfile,
             txt.columns=print.columns,
             txt.cex=print.cex,
             ...)
  } else {
      runsum(object,
             modfile=cmdfile,
             listfile = lstfile,
             ##txt.cex=0.7,
             ...)
  }  
}
