############################################################
#                        cSupport.r                        #
# ---------------------------------------------------------#
# This file contains functions pertaining to the loadC     #
# C compiling/dynamic loading GUI.                         #
#                                                          #
# Authors:                                                 # 
#  Jon T. Schnute <SchnuteJ@pac.dfo-mpo.gc.ca>,            #
#  Anisa Egeli <Anisa.Egeli@dfo-mpo.gc.ca>, and            #
#  Rowan Haigh <HaighR@pac.dfo-mpo.gc.ca>                  #
#                                                          #
############################################################


#loadC----------------------------------2009-07-08
#  Function for launching loadC GUI.
#-----------------------------------------------AE
loadC=function(){
	createWin(system.file("win/loadC.txt", package="PBSmodelling"))
  
  declareGUIoptions("editor")
  readPBSoptions()
  getGUIoptions()
	  
  findPrefix(c("c", "cc", "cpp", "cxx"))
  invisible()
}
#--------------------------------------------loadC


#compileC-------------------------------2009-07-08
#  Compiles a C file into a shared library file, showing and
#  creating a log of the compiler output and possibly alerts
#  if errors occur. If the library is already loaded, it will
#  automatically be unloaded.
# Input:
#  file - the filename of the C file to compile
#  lib - (optional) the name of the resulting shared library
#        file without an extension. If not specified, the
#        prefix of the C file will be used.
#  options - linker options (in one string) to prepend to
#            compilation command}
#  logwindow - If TRUE, the compiler output will be displayed
#     				 in a pop-up window
#  logfile - if TRUE, a log file containing the compiler
#            output will be created
#-----------------------------------------------AE
compileC=function(file, lib="", options="", logWindow=TRUE, logFile=TRUE){
  if(!file.exists(file)){
    showAlert(paste("File ", file, " does not exist in the working directory.",
				sep=""))
    return(invisible())
  }
  
  if(lib=="")
		lib=.stripExt(file)
	lib=.libName(lib)

  try(dyn.unload(lib), silent=TRUE)

  command=paste(R.home(), "/bin/", "R CMD SHLIB -o ", lib, " ", file, " ",
			options, sep="")
  output=system(command, intern=TRUE)
  if(!length(output))
    showAlert(paste(lib, "is already compiled."))
  else{
    output=paste(output, collapse="\n")
    if(logFile)
    	.showLog(paste(command, output, sep="\n\n"), paste(.stripExt(lib), ".log",
					sep=""), noWindow=!logWindow)
		else
			.showLog(paste(command, output, sep="\n\n"), noWindow=!logWindow)
  }
	return(invisible())
}
#-----------------------------------------compileC


#.libName-------------------------------2009-07-08
#  Given a character vector of shared library object
#  names, returns the filenames with the appropriate
#  extension for the user's platform (.dll for Windows or .so
#  for Unix)
# Input:
#  lib - vector of filenames without extensions
# Output:
#   what the corresponding filenames should be on the current
#   platform
#-----------------------------------------------AE
.libName=function(lib=""){
  if (.Platform$OS.type=="windows")
    return(paste(lib, ".dll", sep=""))
  else
    return(paste(lib, ".so", sep=""))
}
#-----------------------------------------.libName


#.guiSource-----------------------------2009-07-08
#  Sources the .r file in the working directory indicated by
#  the prefix entry widget in the GUI
#-----------------------------------------------AE
.guiSource=function(){
  prefix=.getPrefix()
  filename=paste(prefix, ".r", sep="")
  res=try(source(filename))
  if(class(res)=="try-error")
    showAlert(paste("Error sourcing ", filename, ".", sep=""))
}
#---------------------------------------.guiSource


#.guiCompileC---------------------------2009-07-08
#  Gets the prefix and libPrefix arguments from the GUI and
#  uses them to call .compileC
#-----------------------------------------------AE
.guiCompileC=function(){
  prefix=.getPrefix()
  if (is.null(prefix))
    return()

	filename=paste(prefix, ".c", sep="")
	cFiles=Sys.glob(paste(prefix, c("c", "c?", "c??"), sep="."))
	if(length(cFiles) && !any(cFiles==filename))
		filename=cFiles[1]

  getWinVal(c("libPrefix"), scope="L") 
  if (libPrefix==""){
    libPrefix=prefix
    setWinVal(list("libPrefix"=libPrefix))
  }
    
  compileC(filename, libPrefix)
}
#-------------------------------------.guiCompileC


#.guiDyn--------------------------------2009-07-08
#  Based on the previous GUI action, either tries to load or
#  unload the library with the lib prefix specified in the
#  GUI (or the project file prefix if this is left blank
#-----------------------------------------------AE
.guiDyn=function(){
  getWinVal("libPrefix", scope="L")
  if(libPrefix=="")
    libPrefix=.getPrefix()
  if(is.null(libPrefix))
    return()
  
  lib=.libName(libPrefix)
  
  action=getWinAct()[1]
  if(action=="load"){
    if(!file.exists(lib)){
        showAlert(paste("Cannot find", lib, "in working directory"))
        return()
    }
    dyn.load(lib)
  }
  else
    try(dyn.unload(lib), silent=TRUE)
}
#------------------------------------------.guiDyn


#.cleanLoadC----------------------------2009-07-08
#  Clean function
#-----------------------------------------------AE
.cleanLoadC=function(){
	cleanPrefix=.getPrefix(quiet=TRUE)
	if(is.null(cleanPrefix))
		cleanPrefix="*"
	cleanProj(cleanPrefix, suffix=c(".d", ".o", "_res.o", "_res.rc", ".log",
			.libName()), files="Makedeps")
}
#--------------------------------------.cleanLoadC


#.loadCRunComparison--------------------2009-07-08
#  Runs the provided C and R functions a number of times
#  specified in the GUI and writes into text boxes the
#  elapsed time for each
#-----------------------------------------------AE
.loadCRunComparison=function(){
  prefix=.getPrefix()
  if(is.null(prefix))
    return()
    
  rFun=paste(prefix, ".R", sep="")
  cFun=paste(prefix, ".C", sep="")
  if(!exists(rFun) || class(get(rFun))!="function"){
    showAlert(paste("Cannot find function", rFun))
    return()
  }
  if(!exists(cFun) || class(get(cFun))!="function"){
    showAlert(paste("Cannot find function", cFun))
    return()
  }
  getWinVal("runs", winName="loadCGUI", scope="L")
  if (runs<1){
    showAlert("Invalid number of run times.")
    return()
  }
  
  initFun=paste(prefix, ".init", sep="")
  if(exists(initFun) && class(get(initFun))=="function")
    do.call(initFun, list())
  
  gc()
  rTime=proc.time()
  rRet=do.call(rFun, list())
  rTime=proc.time()-rTime
	if(runs>1){
		for(i in 1:(runs-1))
      rTime=rTime+system.time(do.call(rFun, list()), gcFirst=FALSE)
	}

	gc()
  cTime=proc.time()
  cRet=do.call(cFun, list())
  cTime=proc.time()-cTime
  if(runs>1){
    for(i in 1:(runs-1))
      cTime=cTime+system.time(do.call(cFun, list()), gcFirst=FALSE)
  }
  cat("Result for R function:\n")
  print(rRet)
  cat("Result for C function:\n")
  print(cRet)

  setWinVal(list("rTime"=rTime[3], "cTime"=cTime[3], "ratio"=rTime/cTime))
}
#------------------------------.loadCRunComparison


#===== THE END ===================================

