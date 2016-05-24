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

"tabulate.parameters"  <- function(object,prompt=FALSE,outfile=NULL,dir="")
{
  if(prompt==TRUE){
    listfile=paste("run",object@Runno,".lst",sep="")
    ## Get the name of the list file to use
    cat("Type the name of the output file (0=cancel, return=",
        listfile,")\n",sep="")
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
        lstfile <- listfile
      }
    }

  } else {
    lstfile = paste(dir,"run",object@Runno,".lst",sep="")
  }

  if(is.null(lstfile)) {
    cat("The specified file couldn't be found in the current directory.\n")
    return()
  }
  
  parameter.list <- create.parameter.list(lstfile)

  #attach(parameter.list,warn.conflicts=F)  

  ## Set up matrix
    if(any(parameter.list$separval!="" & parameter.list$separval!=0)) {
    ret.mat <- matrix(0,
                      nrow=length(parameter.list$parval),
                      ncol=3,
                      dimnames=list(c(),c("Parameter","Value","RSE"))
                      )
    ret.mat[,1] <- parameter.list$parnam
    ret.mat[,2] <- parameter.list$parval
    ret.mat[,3] <- parameter.list$separval

  } else {
    ret.mat <- matrix(0,
                      nrow=length(parameter.list$parval),
                      ncol=2,
                      dimnames=list(c(),c("Parameter","Value"))
                      )
    ret.mat[,1] <- parameter.list$parnam
    ret.mat[,2] <- parameter.list$parval
  }

  class(ret.mat) <- "char.matrix"

  if(prompt==TRUE){
    cat("Would you like to export the table(s) as a text file? n(y)\n")
    ans <- readline()
  } else {
    if (is.null(outfile)){
      ans = "n"
    } else {
      ans = "y"
    }
  }
  
  if(ans != "y") {
    print.char.matrix(ret.mat,col.names=TRUE)
  }
  else {
    if(prompt==TRUE || is.null(outfile)){
      cat("Please type a filename (excluding the .txt extension):\n"
          )
      ans <- readline()
    } else {
      ans <-  outfile
    }
    print(ret.mat, file = paste(ans, ".txt", sep = ""))
  }
  
  #detach(parameter.list)
  
  return(cat(""))  
  
}
