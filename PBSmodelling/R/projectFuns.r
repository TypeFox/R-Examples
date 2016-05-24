############################################################
#                        projFuns.r                        #
# ---------------------------------------------------------#
# This file contains functions for project management GUI  #
# creation.                                                #
#                                                          #
# Authors:                                                 #
#  Jon T. Schnute <SchnuteJ@pac.dfo-mpo.gc.ca>,            #
#  Anisa Egeli <Anisa.Egeli@dfo-mpo.gc.ca>, and            #
#  Rowan Haigh <HaighR@pac.dfo-mpo.gc.ca>                  #
#                                                          #
############################################################


#cleanProj------------------------------2009-04-21
# Anisa's cleanProj function modified for flexibility.
#--------------------------------------------AE/RH
cleanProj=function(prefix, suffix, files) {
	if (missing(suffix)) suffix = character(0)
	if (missing(files))  files  = character(0)
	rowLen = ceiling(sqrt(max(length(suffix), length(files))))
	if (rowLen == 0) return(invisible(FALSE))
	winDesc = c("window name=cleanWindow title=Clean",
		paste("entry name=cleanPrefix value=\"", prefix, "\" label=Prefix ",
			"mode=character width=12 font=\"bold 9\"", sep = ""),
		"label text=\"\n\nSuffixes to Clean\" font=\"bold 9\"", 
		.makeCleanVec("suff", suffix, rowLen), 
		"label text=\"\n\nFiles to Clean\" font=\"bold 9\"", 
		.makeCleanVec("file", files, rowLen), 
		"grid 1 3 relief=groove padx=4 pady=4", 
		"button function=.selectCleanBoxes action=1 text=\"Select All\" padx=4 pady=4", 
		"button function=.selectCleanBoxes action=0 text=\"Deselect All\" padx=4 pady=4", 
		"button function=doAction text=Clean bg=aliceblue padx=4 pady=4 action=\"PBSmodelling:::.doClean(); closeWin(`cleanWindow`)\"")
	createWin(winDesc, astext = TRUE) 
	invisible(TRUE) }
#----------------------------------------cleanProj


#cleanWD--------------------------------2009-02-24
# Clean all potential garbage files.
#-----------------------------------------------RH
cleanWD=function(files){ # Clean all nuisance files
	rowLen = ceiling(sqrt(length(files)))
	if (rowLen == 0) {
		try(closeWin("cleanWD"),silent=TRUE); return(invisible(FALSE)) }
	winDesc = c("window name=cleanWD title=Clean",
		"label text=\"\n\nFiles to Clean\" font=\"bold 9\"",
		.makeCleanVec("file", files, rowLen),
		"grid 1 3 relief=groove padx=4 pady=4", 
		"button function=.selectCleanBoxes action=1 text=\"Select All\" padx=4 pady=4", 
		"button function=.selectCleanBoxes action=0 text=\"Deselect All\" padx=4 pady=4", 
		"button function=doAction text=Clean bg=aliceblue padx=4 pady=4 action=\"PBSmodelling:::.doCleanWD(); closeWin(`cleanWD`)\"")
	createWin(winDesc, astext = TRUE) 
	invisible(TRUE) }
#------------------------------------------cleanWD


#declareGUIoptions----------------------2012-12-04
#  Used to add options that a GUI uses/loads. The widget
#  names in the window description file should match.
# Input:
#  newOptions - a character vector of option names
#--------------------------------------------AE/RH
declareGUIoptions=function(newOptions){
	.initPBSoptions()
	tget(.PBSmod)
	#.PBSmod$.options$.optionsDeclared<<-.mergeVectors(.PBSmod$.optionsDeclared,newOptions)
	.optionsDeclared=.mergeVectors(.PBSmod$.options$.optionsDeclared,newOptions)
	#packList(".optionsDeclared",".PBSmod$.options")
	.PBSmod$.options$.optionsDeclared <- .optionsDeclared
	tput(.PBSmod)
}
#--------------------------------declareGUIoptions


#findPrefix-----------------------------2012-12-20
# Input:
#  suffix - character vector of suffixes to match to a file.
# Output: character vector of files with matching extensions
# -----------------------------------------------------------
findPrefix=function(suffix, path = "." ) {
	if( length( suffix ) > 1 ) {
		ret <- c()
		for( s in suffix )
			ret <- c( ret, findPrefix( s, path ) )
		return( ret )
	}
	#spat=gsub("\\.","\\\\\\.",suffix)                    # wrong: produces three backslashes
	#spat=gsub("\\.",paste("\\\\","\\.",sep=""),suffix)   # possibly correct but unnecessary
	spat=gsub("\\.","\\\\.",suffix)                       # sufficient
	sfiles=list.files( path, pattern=paste(spat,"$",sep=""),ignore.case=TRUE)
	pref=substring(sfiles,1,nchar(sfiles)-nchar(suffix))
	return(pref)
}
#---------------------------------------findPrefix


#findSuffix-----------------------------2012-12-20
findSuffix=function( prefix, path = "." ) {
	if( length( prefix ) > 1 ) {
		ret <- c()
		for( p in prefix )
			ret <- c( ret, findSuffix( p, path ) )
		return( ret )
	}
	#spat=gsub("\\.","\\\\\\.",prefix)                    # wrong: produces three backslashes
	#spat=gsub("\\.",paste("\\\\","\\.",sep=""),prefix)   # possibly correct but unnecessary
	spat=gsub("\\.","\\\\.",prefix)                       # sufficient
	sfiles=list.files(path,pattern=paste("^", spat,sep=""),ignore.case=TRUE)
	pref=substring(sfiles,nchar(prefix) + 1)
	return(pref)
}
#---------------------------------------findSuffix


#getGUIoptions--------------------------2012-12-04
#  Set used option values as specified by declareGUIoptions in
#  GUI from stored values, possibly trying to load new values
#  from a saved options file
#  It is assumed .PBSmod is initialized
# -------------------------------------------AE/RH
getGUIoptions=function(){
	tget(.PBSmod)
	for(i in .PBSmod$.options$.optionsDeclared){
		option=list()
		option[[i]]=.PBSmod$.options[[i]]
		try(setWinVal(option), silent=TRUE)
	}
}
#------------------------------------getGUIoptions


#getYes---------------------------------2012-12-02
#  A pop-up box prompting the user to choose Yes or No
# Input:
#  message - the message to display in the pop-up
#  title - the title of the pop-up
#  icon - icon to use in message box
# Output: TRUE if Yes was chosen or FALSE if No was chosen
#-----------------------------------------------AE
getYes=function(message, title="Choice", icon="question"){
	answer=as.character(tkmessageBox(message=message, title=title, icon=icon, type="yesno"))
	if(answer=="yes") return(TRUE)
	else return(FALSE)
}
#-------------------------------------------getYes


#openExamples---------------------------2012-12-20
#  Open examples from the examples subdirectory
#  of a given package, making copies into the working
#  directory. If files with the same name already exist than
#  the user is prompted with the choice to overwrite.
# Input:
#  package - package that contains the examples
#            (in example subdir)
#  prefix - prefix of example files
#  suffix - suffixes of example files (character vector)
#-----------------------------------------------AE
openExamples=function(package, prefix, suffix){
  if(missing(package) && missing(prefix) && missing(suffix)){
    fromGUI=TRUE
    action=strsplit(getWinAct()[1], ",")[[1]]
    package=action[1]
    prefix=action[2]
    suffix=action[3:length(action)]
  } else
    fromGUI=FALSE

  filepaths=Sys.glob(paste(system.file("examples", package=package), "/",
      prefix, suffix, sep=""))
  filenames=basename(filepaths)

  for(i in 1:length(filenames)){
    if(!file.exists(filenames[i]) || getYes(paste("Overwrite existing ",filenames[i], "?", sep="")))
      file.copy(from=filepaths[i], to=filenames[i], overwrite=TRUE)
  }

  if(fromGUI)
    try(setWinVal(list("prefix"=prefix)), silent=TRUE)

  for(i in filenames)
    .tryOpen(i)
}
#-------------------------------------openExamples


#promptWriteOptions---------------------2012-12-04
#  Prompts user to save options if a change was made since last load
# Input:
#  fname: name of file in which the options will be saved
# -------------------------------------------AE/RH
promptWriteOptions=function(fname=""){
	.initPBSoptions()
	if(.optionsNotUpdated() && getYes("Set declared options to widget values?"))
		setGUIoptions("*")

	tget(.PBSmod)
	# doesn't save if nothing's changed but we changed the fname - or if we modify .PBSmod directly	
	# if(!is.null(.PBSmod$.options$.optionsChanged) ){
	if(fname=="" && !is.null(.PBSmod$.options$.optionsFile)){
		if (getYes(paste("Save settings to", .PBSmod$.options$.optionsFile, "?")))
			writePBSoptions(.PBSmod$.options$.optionsFile)
	}
	if(fname=="")
		fname="PBSoptions.txt"
	if (getYes(paste("Save settings to", fname, "in working directory?")))
		writePBSoptions(fname)
}
#-------------------------------promptWriteOptions


#setFileOption--------------------------2012-12-20
#  Set a PBS option by browsing for a file.
#  If this is used in a Window description file, the widget
#  action will be used for the option name. If a window
#  containing a widget with the same name as the option is
#  open, the widget value will be set to the new value.
# Input:
#  option - name of the option to change
# Output:
#   Returns TRUE if the option was set, FALSE if the
#   prompt was cancelled
#-----------------------------------------------AE
setFileOption=function(option){
  .setOption(option, "file")
}
#------------------------------------setFileOption


#setwdGUI-------------------------------2012-12-20
#  change the working directory via a GUI
#-----------------------------------------------AE
setwdGUI=function(){
  wd=as.character(tkchooseDirectory())
  if(!length(wd))
    return()
  setwd(wd)
}
#-----------------------------------------setwdGUI


#setGUIoptions--------------------------2012-12-04
#  Transfer option from GUI to option stored in memory.
#  If no option is specified, the option to update is
#  indicated by the last GUI action.
#  "*" updates all options, as given by declareGUIoptions.
# Input:
#  option - name of option to update or "*" to update all
#--------------------------------------------AE/RH
setGUIoptions=function(option){
  .initPBSoptions()
  tget(.PBSmod)
  if(missing(option))
    option=getWinAct()[1]
  if(option=="*"){
    option=names(getWinVal(.PBSmod$.options$.optionsDeclared))
  }
  for(i in option){
    setPBSoptions(i, getWinVal(i)[[i]])
  }
}
#------------------------------------setGUIoptions


#setPathOption--------------------------2012-12-20
#  Set a PBS option by browsing for a directory.
#  If this is used in a Window description file, the widget
#  action will be used for the option name. If a window
#  containing a widget with the same name as the option is
#  open, the widget value will be set to the new value.
# Input:
#  option - name of the option to change
# Output:
#   Returns TRUE if the option was set, FALSE if the
#   prompt was cancelled
#-----------------------------------------------AE
setPathOption=function(option){
  .setOption(option, "dir")
}
#------------------------------------setPathOption


#showAlert------------------------------2012-12-20
#  Show an alert pop-up box with a message.
# Input:
#  message - the message to display in the alert
#  title - the title of the alert box
#  icon - icon to show in alert box
#-----------------------------------------------AE
showAlert=function(message, title="Alert", icon="warning"){
  tkmessageBox(message=message, title=title, icon=icon)
}
#----------------------------------------showAlert


#=================================================
#              HIDDEN FUNCTIONS
#=================================================


#.doClean-------------------------------2009-02-24
# Used by cleanProj(); function called when Clean button is pressed.
#--------------------------------------------AE/RH
.doClean=function(){
	prefix=getWinVal("cleanPrefix",winName="cleanWindow")[[1]]
	vecList=.removeFromList(getWinVal(winName="cleanWindow"), "cleanPrefix")
	filenames=character(0)
	for(i in names(vecList)){
		ii=vecList[[i]] # named logical vector
		type=sub("[[:digit:]]*$", "", i)
		if(type=="suff")
			filenames=c(filenames, Sys.glob(paste(prefix,names(ii)[ii], sep="")))
		else
			filenames=c(filenames, Sys.glob(names(ii)[ii]))
#if (i=="suff3") {browser();return()}
	}
	if(!length(filenames))
		showAlert("No files to delete.")
	else if(getYes(paste("Delete ", paste(filenames, collapse=", "), "?",sep="")))
		file.remove(filenames) 
	remaining=file.exists(filenames)
	if(sum(remaining)) 
		showAlert(paste("Failed to delete",paste(filenames[remaining],collapse=", ")))
}
#-----------------------------------------.doClean


#.doCleanWD-----------------------------2009-02-24
# Anisa's .doClean function modified for file names only
#-----------------------------------------------RH
.doCleanWD=function () { 
	vec=getWinVal(winName="cleanWD",scope="L")
	vecList=logical()
	for (i in names(vec)) vecList=c(vecList,vec[[i]])
	filenames = names(vecList)[vecList]
	filenames=Sys.glob(filenames)
	if (!length(filenames)) 
		showAlert("No files to delete.")
	else if (getYes(paste("Delete ", paste(filenames, collapse = ", "), "?", sep = ""))) 
		file.remove(filenames)
	remaining = file.exists(filenames)
	if (sum(remaining)) 
		showAlert(paste("Failed to delete", paste(filenames[remaining], collapse = ", "))) 
}
#---------------------------------------.doCleanWD


#One big hack to allow GUIs from outside PBSmodelling to call PBSmodelling's internal hidden functions.
#Do not rely on this to remain here - functions beginning with a '.' should ONLY be called from PBSmodelling.
#TODO REMOVE THIS - after PBSadmb is refactored
.getHiddenEnv <- function()
{
	return( environment() )
}


#.getHome-------------------------------2012-12-20
#  Returns platform dependent home drive (Windows) or user
#  home (Unix) -- NOT TESTED ON UNIX.
# Output:
#  HOMEDRIVE or HOME environment variable for Windows or Unix
#  respectively
#-----------------------------------------------AE
.getHome=function(){
  if (.Platform$OS.type=="windows")
    return(Sys.getenv("HOMEDRIVE")[[1]])
   return(Sys.getenv("HOME")[[1]])
}
#-----------------------------------------.getHome


#.getPrefix-----------------------------2012-12-20
#  get the project prefix from the focused GUI. Used to
#  standardize the error message.
# Input:
#  quiet - If TRUE, no errors/alerts will be displayed
# Output: the prefix, or NULL if the entry box was empty.
#-----------------------------------------------AE
.getPrefix=function(quiet=FALSE){
  getWinVal("prefix", scope="L")
  if (prefix==""){
    if(!quiet)
      showAlert("Please enter the prefix/name of your project.")
    return(NULL)
  }
  return(prefix)
}
#---------------------------------------.getPrefix


#.mergeLists----------------------------2012-12-20
#  Add a second list to a first list. If both lists share any
#  keys, the values from the second list overwrite those in
#  the first list. Order of components in list is preserved.
# Input:
#  list1 - first list
#  list2 - second list
# Output:
#   the two lists merged as described above
#-----------------------------------------------AE
.mergeLists=function(list1, list2){
  if(!length(list1))
    return(list2)
  if(!length(list2))
    return(list1)

  newComponents=list()

  for(i in 1:length(list2)){
    nameI=names(list2[i])
    if(any(names(list1)==nameI))
      list1[[nameI]]=list2[[nameI]]
    else
      newComponents[[nameI]]=list2[[nameI]]
  }
  return(c(list1, newComponents))
}
#--------------------------------------.mergeLists


#.mergeVectors--------------------------2012-12-20
#  Add a second vector to a first vector. If both vectors
#  share any values, the resulting vector will only contain
#  this value once
# Input:
#  v1 - first vector
#  v2 - second vector
# Output:
#   the two vectors merged as described above
#-----------------------------------------------AE
.mergeVectors=function(v1, v2){
  if(!length(v2))
    return(v1)

  newVals=vector()
  for(i in 1:length(v2)){
    if(!any(v1==v2[i]))
      newVals=c(newVals, v2[i])
  }
  return(c(v1, newVals))
}
#------------------------------------.mergeVectors


#.optionsNotUpdated---------------------2012-12-04
#  Check if any of the options given by declareGUIoptions are
#  different in the GUI from their stored values. A blank
#  GUI entry is equivalent to the option not being set.
# Ouput:
#  Returns TRUE if any of the options are different in the
#  GUI than in their stored values
# -------------------------------------------AE/RH
.optionsNotUpdated=function(){
	.initPBSoptions()
	tget(.PBSmod)
	if(is.null(.PBSmod$.options$.optionsDeclared))
		return(FALSE)
	winVals=try(getWinVal(.PBSmod$.options$.optionsDeclared), silent=TRUE)
	if(class(winVals)=="try-error")
		return(FALSE)
	for(i in names(winVals)){
		if((is.null(.PBSmod$.options[[i]]) && winVals[[i]]!="") || 
		(!is.null(.PBSmod$.options[[i]]) && winVals[[i]]!=.PBSmod$.options[[i]]))
			return(TRUE)
	}
	return(FALSE)
}
#-------------------------------.optionsNotUpdated


#.removeFromList------------------------2012-12-20
#  Remove components of a list with the given names. "NA" can
#  be used to remove NA names
# Input:
#  l - a list
#  items - character vector of names of components to remove
# Output:
#  a list with the chosen components removed
#-----------------------------------------------AE
.removeFromList=function(l, items){
  if(!length(l) || !length(items))
    return(l)

  indices=integer()
  for(i in 1:length(l)){
    nameI=names(l[i])
    if(is.na(nameI))
      nameI="NA"
    if(any(items==nameI)){
      indices=c(indices, i)
    }
  }
  if(length(indices))
    return(l[-indices])
  return(l)
}
#----------------------------------.removeFromList


#.selectCleanBoxes----------------------2008-08-25
# This is used by cleanProj. It is the function for selecting or deselecting
# all of the checkboxes.
#-----------------------------------------------AE
.selectCleanBoxes=function(){
	action=getWinAct()[1]
	vecList=.removeFromList(getWinVal(), "cleanPrefix")
	for(i in 1:length(vecList)){
		vecList[[i]]=rep(action, length(vecList[[i]]))
	}
	setWinVal(vecList)
}
#--------------------------------.selectCleanBoxes


#.setOption-----------------------------2012-12-04
# This used used for setPathOption and setFileOption
# See these functions for a description. 
#The type is "dir" or "file" for each of those functions.
#--------------------------------------------AE/RH
.setOption=function(option, type){
	if(missing(option))
		option=getWinAct()[1]
	.initPBSoptions()
	tget(.PBSmod)
	if(!is.null(.PBSmod$.options[[option]]) && file.exists(.PBSmod$.options[[option]]))
		initDir=.PBSmod$.options[[option]]
	else
		initDir=.getHome()
	if(type=="dir")
		value=paste(as.character(tkchooseDirectory(initialdir=initDir)), collapse=" ")
	else
		value=paste(as.character(tkgetOpenFile(initialdir=initDir)), collapse=" ")
	if(nchar(value)){
		newVal=list()
		newVal[[option]]=value
		try(setWinVal(newVal), silent=TRUE)
		setPBSoptions(option, value)
		return(TRUE)
	}
	return(FALSE)
}
#---------------------------------------.setOption


#.showLog-------------------------------2012-12-20
#  Given output for a log, will make pop-up log window
#  containing this output and/or write this output to a
#  logfile. \r is stripped from the log text for display
# Input:
#  logText - the text for the log
#  fname (optional) - name for log file
#  noWindow - if TRUE, log window will not be shown
#  width - width of log window
#  height - height of log window
#-----------------------------------------------AE
.showLog=function(logText, fname, noWindow=FALSE, width=80, height=30){
  if(!noWindow){
    winDesc=c(
        "window name=PBSlog title=LOG",
        paste("text name=logText width=", width, " height=", height,
        " edit=FALSE scrollbar=TRUE mode=character", sep=""))
    createWin(winDesc, astext=TRUE)
    setWinVal(list("logText"=gsub("\r", "", logText)), winName="PBSlog")
  }
  if(!missing(fname)){
    cat(logText, file=fname)
  }
}
#-----------------------------------------.showLog


#.stripExt------------------------------2012-12-20
#  remove file extension from end of filename
# Input:
#  x - character vector of filenames, ie. "foo.c"
# Output: the filename without the extension, ie. "foo"
#-----------------------------------------------AE
.stripExt=function(x){
	return(sub("[.].{1,3}$", "", x))
}
#----------------------------------------.stripExt


#.tryOpen-------------------------------2012-12-04
#  Tries to open a given file using an editor entered by the
#  GUI. If an editor wasn't set, tries to open using openFile.
#  Appropriate alerts are shown if quiet isn't turned on.
# Input:
#  filename - the file to open
#  quiet - if quiet is TRUE, alerts will not be shown
# Output:
#   returns TRUE if file exists and a program to open the
#   file was specified somewhere, FALSE otherwise.
# -------------------------------------------AE/RH
.tryOpen=function(filename, quiet=FALSE){
	filename=filename[1]
	if(!is.character(filename) || !file.exists(filename)){
		if(!quiet)
			showAlert(paste("File", filename, "does not exist."))
		return(FALSE)
	}
	.initPBSoptions()
	tget(.PBSmod)
	if(!is.null(.PBSmod$.options$editor) && file.exists(.PBSmod$.options$editor)){
		cmd=paste(shQuote(.PBSmod$.options$editor), filename)
		system(cmd, wait=FALSE, invisible=FALSE)
	} else {
		tryRet=try(openFile(filename), silent=TRUE)
		if(class(tryRet)=="try-error"){
			if(!quiet)
				showAlert(paste("Could not open file ", filename, ". Please choose ",
				"a default editor or set the appropriate file association.", sep=""))
		return(FALSE)
		}
	}
	return(TRUE)
}
#-----------------------------------------.tryOpen


#.makeCleanVec--------------------------2009-03-03
# This is used by cleanProj() to create the strings describing checkbox vectors.
#--------------------------------------------AE/RH
.makeCleanVec=function(vecName, items, rowLen){
	vecDesc=character(0)
	nItems=length(items)
	nVecs=ceiling(nItems/rowLen)
	if(!nVecs)
		return(vecDesc)
	for(i in 1:nVecs){
		vecI=paste("vector names=", vecName, i, " ", sep="")
		currItems=items[((i-1)*rowLen+1):min(rowLen*i, nItems)]
		itemStr=paste("\'",paste(currItems, collapse="\' \'"),"\'",sep="")
		vecI=paste(vecI, 'vecnames="', itemStr, '" ', sep="")
		vecI=paste(vecI, 'labels="', itemStr, '"', sep="")
		vecI=paste(vecI, " length=", length(currItems), sep="")
		vecI=paste(vecI," mode=logical vertical=FALSE value=TRUE padx=4 pady=4")
		vecDesc=c(vecDesc, vecI)
	}
	return(vecDesc)
}
#------------------------------------.makeCleanVec


#===== THE END ===================================

