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
gofSetup <- function(runno,save,onefile,saveType,pageWidth,pageHeight) {

  ## Determine if we are running a Mac so that we can behave appropriately.
  isWin <- FALSE
  isMac <- FALSE
  isLin <- FALSE
  if(length(grep("darwin",R.version$os))!=0) { isMac <- TRUE}
  if(length(grep("linux",R.version$os))!=0)  { isLin <- TRUE}
  if(length(grep("mingw",R.version$os))!=0)  { isWin <- TRUE}
  
  ## Get access to libraries that we need
  #library(lattice)
  #library(grid)
  #library(Hmisc)
  #library(xpose4)

  ## Create file name to store graphs in. The graphs wiíll be saved in
  ## a sub-directory called Plots, which will be created if it doesn't exist
  saveDir <- "Plots/"
  if(is.na(file.info(saveDir)$isdir)) dir.create(saveDir)

  ## Create the file name to store graphs in
  runnoLabel   <- ""
  if(!is.null(runno)) runnoLabel <- paste("run",runno,sep="")
  saveFileName                   <- paste(saveDir,sys.calls()[[1]],runnoLabel,
                                          "-%02d.",saveType,sep="")

  ## Start a device, with recording, the same size as the save file size if we are not saving
  if(!save) {
    if(isMac) trellis.device(device="quartz",  new = FALSE,width=pageWidth,height=pageHeight)
    if(isWin) trellis.device(device="windows", new = FALSE,width=pageWidth,height=pageHeight,
                             record=TRUE)
    if(isLin) trellis.device(new = FALSE,width=pageWidth,height=pageHeight)
  } else {
    if(saveType=="png") png(filename=saveFileName,width=pageWidth,height=pageHeight,res=150,
         units="in",bg="white")
    if(saveType=="pdf") pdf(file=saveFileName,width=pageWidth,height=pageHeight,
         onefile=onefile,bg="white")
    if(saveType=="wmf") eval(parse(text=paste("win.metafile(filename=",saveFileName,
                                       ",width=",pageWidth,",height=",pageHeight,")",sep="")))
  }
}
