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

read.vpctab <- function(vpctab=NULL,
                        object=NULL,
                        vpc.name="vpctab",
                        vpc.suffix="",
                        tab.suffix="",
                        inclZeroWRES=FALSE,
                        verbose=FALSE,
                        ...) {
    ## Make sure we have the necessary variables defined
    if(is.null(object) & is.null(vpctab)){
        cat(paste("Both the arguments object and vpctab are NULL\n"))
        cat(paste("At least one of these must be defined\n"))
        return(NULL)
    }
    if(is.null(vpctab)){
        vpc.file <- sapply(vpc.name,paste,object@Runno,vpc.suffix,tab.suffix,sep="")
    } else {
        vpc.file <- vpctab
    }
    
    if(!is.readable.file(vpc.file)) {
        if (verbose) cat(paste(vpc.file,"not readable\n"))
        return(NULL)
    } else { 
        if(verbose) cat(paste("    Reading",vpc.file,"\n"))
        orig.data <- read.table(file=vpc.file,header=T,sep=",",colClasses="numeric")
    }
    
    ## check that classes are present
    if (length(findClass("xpose.data")) < 1) {
      createXposeClasses()
    }

    ## get the VPC run number
    if(!is.null(vpctab)){
      vpctabnum.start <- gregexpr("vpctab[[:digit:]]",vpctab)
      if(vpctabnum.start!=-1){
        vpc.num <- substring(vpctab,vpctabnum.start[[1]][1]+6)
      } else {
        vpc.num <- "0"
      }
    } else {
      vpc.num <- object@Runno
    }
          
    ##create vpcdb
    vpcdb       <- new("xpose.data",
                       Runno=vpc.num,
                       Doc=NULL,
                       Data = NULL 
                       )

    ## read local options
    if (is.readable.file("xpose.ini")) {
      vpcdb <- xpose.read(vpcdb, file="xpose.ini")
    } else {
      ## read global options
      rhome   <- R.home()
      xdefini <- paste(rhome, "\\library\\xpose4\\xpose.ini", sep="")
      if (is.readable.file(xdefini)) {
        vpcdb <- xpose.read(vpcdb, file=xdefini)
      }
    }
    
    ## Add data to vpcdb
    Data(vpcdb) <- orig.data

    
    ## if object exists add vpcdb to @Data and @Labels structure
    if(!is.null(object)){
      if(!inclZeroWRES) {
        zero.wres.locs <- object@Data[,xvardef("wres",object)]!=0
      } else {
        zero.wres.locs <- rep(TRUE,length(object@Data$ID))
      }
      obj.data.red  <- Data(object,inclZeroWRES,onlyfirst=FALSE,subset=NULL)

      ## check that the data objects have the same length
      if(!length(obj.data.red$ID) == length(vpcdb@Data$ID)){
        cat(paste(vpc.file,"and the data in the xpose database ('object')\n"))
        cat(paste("have different lengths after subsetting on WRES\n"))
        cat(paste("This must be resolved\n"))
        return(NULL)
      }
      ## check that the data objects have the same ID values
      if(!all(obj.data.red$ID == vpcdb@Data$ID)){
        cat(paste(vpc.file, "and the data in the xpose database ('object')\n"))
        cat(paste("have different ID values on some rows in the data after subsetting on WRES\n"))
        cat(paste("This must be resolved\n"))
        return(NULL)
      }

      ## merge the data
      vpc.names <- names(vpcdb@Data)
      obj.names <- names(object@Data)
      extra.names <- c()
      for(i in vpc.names){
        tmp <- grep(i,obj.names)
        if((length(tmp)==0)) extra.names <- c(extra.names,i)
      }
      object@Data[extra.names] <- NA
      object@Data[zero.wres.locs,vpc.names] <- vpcdb@Data
      
      ## add extra names to Labels list
      object@Prefs@Labels[extra.names] <- extra.names
    } else { 
      object <- vpcdb
      ##inclZeroWRES=TRUE
    }

    return(object)

  }
