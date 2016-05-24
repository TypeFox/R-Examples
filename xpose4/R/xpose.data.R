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

xpose.data <-function(runno,
                      tab.suffix="",
                      sim.suffix="sim",
                      cwres.suffix="",
                      directory="",
                      quiet=TRUE,
                      table.names=c("sdtab","mutab","patab","catab",
                        "cotab","mytab","extra","xptab","cwtab"),
                      cwres.name=c("cwtab"),
                      mod.prefix="run",
                      mod.suffix=".mod",
                      phi.suffix=".phi",
                      phi.file=NULL,
                      ##vpc.name="vpctab",
                      nm7=NULL,  # T/F if table files are for NM7, NULL for undefined
                      ...) {


  ##options(warn=-1) # suppress warnings

  ## make table lists
  match.pos <- match(cwres.name,table.names)
  if (!is.na(match.pos)) table.names <- table.names[-match.pos]

  ## Create the table file names to process
  myfun <- function(x,directory,runno,cwres.suffix,sim.suffix,tab.suffix) {
    paste(directory,x,runno,cwres.suffix,sim.suffix,tab.suffix,sep="")
  }

  tab.files   <- sapply(table.names,myfun,directory,runno,cwres.suffix="",sim.suffix="",tab.suffix)

  cwres.files <- sapply(cwres.name,myfun,directory,runno,cwres.suffix,sim.suffix="",tab.suffix)

  sim.files   <- sapply(table.names,myfun,directory,runno,cwres.suffix="",sim.suffix,tab.suffix)

  cwres.sim.files <- sapply(cwres.name,myfun,directory,runno,cwres.suffix,sim.suffix,tab.suffix)

  tab.files <- c(tab.files,cwres.files)
  sim.files <- c(sim.files,cwres.sim.files)

  ## Read the table files.
  cat("\nLooking for NONMEM table files.\n")
  tmp <- read.nm.tables(table.files=tab.files,
                        quiet=quiet,...)

  ## Fail if we can't find any.
  if(is.null(tmp)) {
    cat("Table files not read!\n")
    return(NULL)
  }

  ## check if NM.version is > 6
  if(is.null(nm7)){
    if(any(!is.na(match(c("IPRED","IWRES"),names(tmp))))){
      nm7 <- T
    } else {
      nm7 <- F
    }
  }


  ## check that classes are present
  #if (length(findClass("xpose.data")) < 1) {
  #  createXposeClasses(nm7=nm7)
  #}
  createXposeClasses(nm7=nm7)

  ## Create the object
  xpobj       <- new("xpose.data",
                     Runno=runno,
                     Doc=NULL,
                     Data = NULL #read.nm.tables(runno,tab.suffix=tab.suffix,
                       #quiet=TRUE)
                     )

  ## read local options
  if (is.readable.file("xpose.ini")) {
    xpobj <- xpose.read(xpobj, file="xpose.ini")
  } else {
    ## read global options
    rhome   <- R.home()
    xdefini <- paste(rhome, "\\library\\xpose4\\xpose.ini", sep="")
    if (is.readable.file(xdefini)) {
      xpobj <- xpose.read(xpobj, file=xdefini)
    }else{
      xdefini2 <- paste(rhome, "\\library\\xpose4\\xpose.ini", sep="")
      if (is.readable.file(xdefini2)) {
        xpobj <- xpose.read(xpobj, file=xdefini2)
      } 
    }
  }

  ## read tmp data into xpose object
  Data(xpobj) <- tmp
  cat("Table files read.\n")


  ## read phi file
  ind.data <- NULL
  nsim.phi <- NULL
  if(nm7){
    phi.data <- read.phi(phi.file=phi.file,
                         phi.prefix=mod.prefix,
                         runno=runno,
                         phi.suffix=phi.suffix,
                         ##sim.suffix="sim",
                         quiet=quiet,
                         nm7=nm7,
                         directory=directory,
                         ...)
    # browser() Elins bug
    if(!is.null(phi.data)){
      ## check that phi file has right size
      if(dim(phi.data)[1]==dim(unique(xpobj@Data[xvardef("id",xpobj)]))[1]){
        xpobj@Data.firstonly <- phi.data
      }else{
        ## get the first unique ID values data from phi file
        first.phi.data <- phi.data[!duplicated(phi.data[,xvardef("id",xpobj)]),]
        sim.phi.data <- phi.data[duplicated(phi.data[,xvardef("id",xpobj)]),]

        xpobj@Data.firstonly <- first.phi.data

        nsim.phi.nrows <- dim(sim.phi.data)[1]
        first.phi.nrows <- dim(first.phi.data)[1]
        if(regexpr("\\.",as.character(nsim.phi.nrows/first.phi.nrows)) !=-1) {
          cat("The length of the Phi data and the Phi simulated data do not match!\n")
          return(xpobj)
        }
        nsim.phi <- nsim.phi.nrows/first.phi.nrows
      }
    }
  }
  

  
  




  ## is.readable loop
  ## check orig files + sim is.readable output & sim must match


  ## error handling for simulations!
  cat("\nLooking for NONMEM simulation table files.\n")
  gosim <- TRUE
  simct <- FALSE

  ## check if there are any simulation files
  for(i in 1:length(sim.files)) {
      if (is.readable.file(sim.files[i]))  {
          simct <- TRUE
      }
  }

  if (simct){
    for(i in 1:length(tab.files)) {
      if ((is.readable.file(tab.files[i])) && (!is.readable.file(sim.files[i])))  {
        err.mess <- paste(sim.files[i],"not found!")
        gosim <- FALSE
        break
      }
    }
  } else {
    gosim <- FALSE
  }


  if (gosim==FALSE) {

    if (!simct) {
      #cat("  Files are either not present or not named correctly\n")
      #cat("  (e.g. sdtab1a instead of sdtab1sim)\n")
    } else {
      cat("  There is not the same number of normal and \n")
      cat("  simulation table files for the current run number:\n")
      cat(paste("  ",err.mess,"\n",sep=""))
    }
    cat("No simulated table files read.\n\n")
  }

  if (gosim==TRUE) {
    simtmp <- read.nm.tables(sim.files,
                             #runno,
                             #tab.suffix=paste(sim.suffix,tab.suffix,sep=""),
                             #cwres.suffix=paste(sim.suffix,cwres.suffix,sep=""),
                             quiet=quiet)

    if(!is.null(tmp)) {
      SData(xpobj) <- simtmp
      cat("Simulation table files read.\n")
    } else {
      cat("There was a problem reading the simulation tables!\n")
      cat("Simulation tables not read!\n")
      return(NULL)
    }

    if(!is.null(nsim.phi)){
      ## check that there are the same number of simulations in the phi file and the table files
      if(!(xpobj@Nsim == nsim.phi)){
        cat("\nThere are not the same number of simulations\n",
            "in the table files and the phi file.\n",
            "Something is wrong with the phi file.\n",
            "It will not be used.\n",sep="")
        xpobj@Data.firstonly <- NULL
      } else {
        xpobj@SData.firstonly <- sim.phi.data
      }
    }
  }


  ## read local options
  if (is.readable.file("xpose.ini")) {
    xpobj <- xpose.read(xpobj, file="xpose.ini")
  } else {
    ## read global options
    rhome   <- R.home()
    xdefini <- paste(rhome, "\\library\\xpose4\\xpose.ini", sep="")
    if (is.readable.file(xdefini)) {
      xpobj <- xpose.read(xpobj, file=xdefini)
    }
  }

  ## clean up
  if(file.exists(".sdtab.names.tmp")) file.remove(".sdtab.names.tmp")
  if(file.exists(".catab.names.tmp")) file.remove(".catab.names.tmp")
  if(file.exists(".cotab.names.tmp")) file.remove(".cotab.names.tmp")
  if(file.exists(".patab.names.tmp")) file.remove(".patab.names.tmp")

  tmp.obj <- read.vpctab(object=xpobj,
                         #inclZeroWRES=inclZeroWRES,
                         tab.suffix=tab.suffix,
                         ...)
  if(!is.null(tmp.obj)) xpobj <- tmp.obj

  if(is.null(check.vars(c("idv"),xpobj))) {
    cat("\n*********PLEASE NOTE: idv NOT IDENTIFIED************\n")
    cat("The independent variable (idv) has not been identified\n")
    cat("in the table files!  Please use the command line function\n")
    cat("'change.xvardef' (use '?change.xvardef' for help) or the classic\n")
    cat("menu system (select: Preferences/Manage variables/Change idv)\n")
    cat("to identify the name of the idv\n")
    cat("****************************************************\n")
  }



  return(xpobj)


}

