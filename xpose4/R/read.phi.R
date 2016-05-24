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

read.phi <-
  function(phi.file=NULL,
           phi.prefix="run",
           runno=NULL,
           phi.suffix=".phi",
           ##sim.suffix="sim",
           quiet=TRUE,
           nm7=TRUE,
           directory="",
           ...)
{
  
  if (!nm7){
    if(!quiet) cat("This function only works for NONMEM 7")
    return(NULL)
  }

  ## create file name
  if(is.null(phi.file)){
    if(is.null(runno)) {
      cat(paste("runno must be specified if no phi file name is provided\n"))
      return(NULL)
    }
    filename <- paste(directory,phi.prefix,runno,phi.suffix,sep="")
  } else {
    filename <- phi.file
  }
  
  if(!is.readable.file(filename)) {
    if (!quiet) {cat(filename,"not readable\n")}
    return(NULL)
  } else {
    cat(paste("    Reading",filename,"\n"))
    
    ##ind.vals <- read.table(filename,header=T,skip=1)
    
    ## Check which type of separator we have in our tables
    header.line = scan(file=filename,nlines=1,skip=1,what="character",sep="\n",quiet=T)
    sep.char = ""
    ##if(length(grep(",",header.line))!=0) sep.char = ","
    
    ## Check if we have unequal number of fields in the file
    ## used for multiple simulations
    fields.per.line <- count.fields(filename)
    fields.in.first.line <- fields.per.line[1]
    fields.in.rest <- fields.per.line[-1]
    if((length(unique(fields.in.rest))!=1) ||
       (all(fields.in.first.line==fields.in.rest))){ 
      if(!quiet) {
        cat(paste("Found different number of columns in different rows of ",filename,".\n",sep=""))
        cat("This may be due to multiple header rows \n")
        cat("caused by running multiple simulations in NONMEM (NSIM > 1).\n")
        cat("Will try to remove these rows. It may take a while...\n")
      }
      tmp   <- readLines(filename, n = -1)
      inds  <- grep("TABLE",tmp)
      if (length(inds)!=1){
        inds  <- inds[c(2:length(inds))]
        inds2 <- inds+1
        tempfile<- paste(filename,".xptmp",sep="")
        write.table(tmp[-c(inds,inds2)],file=tempfile,
                    row.names=FALSE,quote=FALSE)
        ind.vals <- read.table(tempfile,skip=2,header=T,sep=sep.char)
        unlink(tempfile)
      } else {
        ind.vals <- read.table(filename,skip=1,header=T,sep=sep.char)
      }
    } else {
      ind.vals <- read.table(filename,skip=1,header=T,sep=sep.char)
    }
    return(ind.vals)
  }
}
