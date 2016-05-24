## Copyright (c) 2004-2012, Ph. Grosjean <phgrosjean@sciviews.org>
##
## This file is part of ZooImage
##
## ZooImage is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## ZooImage is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with ZooImage.  If not, see <http://www.gnu.org/licenses/>

ZIDlg <- function ()
{
	## In this version, we use a simpler implementation, using svDialogs
	## and menus added to RGui, JGR or ctxmenu
	ZIname <- getTemp("ZIname")
	menuDel(ZIname)
	menuAdd(ZIname)
	menuAddItem(ZIname, "Load objects", "loadObjects()")
	menuAddItem(ZIname, "Save objects", "saveObjects()")
	menuAddItem(ZIname, "List objects", "listObjects()")
	menuAddItem(ZIname, "Remove objects", "removeObjects()")
	menuAddItem(ZIname, "-", "")
	menuAddItem(ZIname, "Online help", 'help("zooimage")')
	menuAddItem(ZIname, "Manual", "viewManual()")
	menuAddItem(ZIname,
		"Web site", 'browseURL("http://www.sciviews.org/zooimage")')
	menuAddItem(ZIname, "--", "")
	menuAddItem(ZIname, "About...", "aboutZI(TRUE)")
	
	menuDel("Analyze")
	menuAdd("Analyze")
	menuAddItem("Analyze", "Acquire images...", "acquireImg()")
	menuAddItem("Analyze", "Import images...", "importImg()")
	menuAddItem("Analyze", "Process images...", "processImg()")
	menuAddItem("Analyze", "Make .zid files...", "makeZid()")
	menuAddItem("Analyze", "-", "")
	menuAddItem("Analyze", "Make training set...", "makeTrain()")
	menuAddItem("Analyze", "Add vignettes to training set", "addVigsToTrain()")
	menuAddItem("Analyze", "Read training set...", "collectTrain()")
	menuAddItem("Analyze", "Make classifier...", "makeClass()")
	menuAddItem("Analyze", "Analyze classifier...", "analyzeClass()")
	menuAddItem("Analyze", "Automatic classification of vignettes...",
		"vignettesClass()") 
	menuAddItem("Analyze", "--", "")
	menuAddItem("Analyze", "Edit samples description", "editDescription()")
	menuAddItem("Analyze", "Process samples...", "processSamples()")
	menuAddItem("Analyze", "View results...", "viewResults()")
	menuAddItem("Analyze", "Export results...", "exportResults()")
	
	## Menu 'Functions' not added yet!
	
	menuDel("Utilities")
	menuAdd("Utilities")
	menuAddItem("Utilities", "Calibrate grayscale (16bit)", "calib()")
	menuAddItem("Utilities", "Biomass conversion specification",
		"fileEdit(file.path(getTemp('ZIetc'), 'Conversion.txt'))")
	menuAddItem("Utilities", "-", "")
	menuAddItem("Utilities", "Image viewer( XnView)", 'startPgm("ImageViewer")')
	menuAddItem("Utilities", "Image analyzer (ImageJ)",
		'startPgm("ImageEditor", switchdir = TRUE, iconize = TRUE)')
	menuAddItem("Utilities", "Metadata editor",
		'fileEdit(selectFile("ZimZis"))')
	menuAddItem("Utilities", "Simple acquisition (Vuescan)",
		'startPgm("VueScan", switchdir = TRUE)')
	menuAddItem("Utilities", "--", "")
	menuAddItem("Utilities", "New R graph", "dev.new()")
	menuAddItem("Utilities", "Activate next graph",
		"{dev.set(); if (isRgui()) bringToTop()}")
	menuAddItem("Utilities", "Close all graphs", "graphics.off()")
	menuAdd("Utilities/Options")
	menuAddItem("Utilities/Options", "Change active dir...",
		"setwd(dlgDir()$res)")
	menuAddItem("Utilities/Options", "-", "")
	menuAddItem("Utilities/Options", "Define decimal separator",
		"optInOutDecimalSep()")
	
#	## This is the old implementation usig svWidgets
#	# If the window is already created, just activate it...
#	if ("ZIDlgWin" %in% WinNames()) {
#		ZIDlgWin <- WinGet("ZIDlgWin")
#		tkfocus(ZIDlgWin)  	# Doesn't work with Rgui.exe, but next command does
#		tkwm.deiconify(ZIDlgWin)
#    	return(invisible())
#	}
#
#	# Construct the window
#	tkWinAdd("ZIDlgWin", title = paste(getTemp("ZIname"), "assistant"),
#		pos = "-100+10")
#	ZIDlgWin <- WinGet("ZIDlgWin")
#
#	# Do not show it until it is completelly constructed!
#	tkwm.withdraw(ZIDlgWin)
#	on.exit(tkwm.deiconify(ZIDlgWin))
#
#	# Change the icon of that window (if under Windows)
#	if (isWin()) tk2ico.set(ZIDlgWin, getTemp("ZIico"))
#
#	# Add a menu (load it from a spec file)
#	Pkg <- getTemp("ZIguiPackage", default = "zooimage")
#	MenuReadPackage(Pkg, file = "MenusZIDlgWin.txt")
#
#	# Add a toolbar (read it from file 'ToolbarZIDlgWin.txt')
#	ToolRead(file.path(getTemp("ZIgui"), "ToolbarsZIDlgWin.txt"))
#
#	# Add a statusbar with a text and a progressbar
#	status <- tk2frame(ZIDlgWin)
#	statusText <- tk2label(status, text = paste("Ready -", getwd()),
#		justify = "left", anchor = "w", width = 60)
#	statusProg <- tk2progress(status, orient = "horizontal", maximum = 100)
#	tkpack(statusProg, side = "right")
#	tkpack(statusText, side = "left", fill= "x")
#	tkpack(status, side = "bottom", fill = "x")
#	tkpack(tk2separator(ZIDlgWin), side = "bottom", fill = "x")
#
#	# Keep track of statusText / statusProg
#	assignTemp("statusText", statusText)
#	assignTemp("statusProg", statusProg)
#	## Change value of the progressbar
#	#tkconfigure(getTemp("statusProg"), value = 50)
#	## Change text of the statusbar
#	#tkconfigure(getTemp("statusText"), text = paste("Ready -", getwd()))

#	## Add a function for progress() to update the progressbar in this window
#	assignTemp(".progress", list(progressZIGUI = function (value, max.value) {
#		if (!"ZIDlgWin" %in% WinNames()) return()
#	
#		if (is.null(max.value)) {
#		    max.value <- 100
#		    percent <- TRUE
#		} else percent <- FALSE
#	
#		if (value > max.value) { # Erase progressmeter
#			rmTemp("statusBusy")
#			tkconfigure(getTemp("statusProg") , value = 0)
#			tkconfigure(getTemp("statusText") , text = paste("Ready -", getwd()))		
#		} else { # Show progress
#			assignTemp("statusBusy", TRUE)
#			## Calculate fraction and show it in the progress bar
#			if (!percent) value <- value / max.value * 100
#			tkconfigure(getTemp("statusProg"), value = value)
#			## Display the progress text also in the statusbar
#			tkconfigure(getTemp("statusText"), text = message)
#		}
#		.Tcl("update idletasks")
#	}))
#
#	if (!isWin()) {
#		# The activate R console & R graph do not work elsewhere
#        MenuStateItem("$Tk.ZIDlgWin/Apps", "&R Console", FALSE)
#		MenuStateItem("$Tk.ZIDlgWin/Apps", "Active R &Graph", FALSE)
#	}
#
#	# For each of the six external programs, look if they are accessible,
#	# otherwise, inactivate
#	if (is.null(getOption("fileEditor")))
#         MenuStateItem("$Tk.ZIDlgWin/Apps", "&Metadata editor", FALSE)
#    if (is.null(getOption("ImageEditor")))
#         MenuStateItem("$Tk.ZIDlgWin/Apps", "Image &analyzer (ImageJ)", FALSE)
#    if (is.null(getOption("ImageViewer")))
#         MenuStateItem("$Tk.ZIDlgWin/Apps", "Image &viewer (XnView)", FALSE)
#    if (is.null(getOption("VueScan")))
#         MenuStateItem("$Tk.ZIDlgWin/Apps", "Simple acquisition (&VueScan)", FALSE)
#
#	# Change the window to non resizable and topmost (f under Windows)
#	if (isWin()) tcl("wm", "attributes", ZIDlgWin, topmost = 1)
#	tkwm.resizable(ZIDlgWin, 0, 0)
#	# Focus on that window
#	tkfocus(ZIDlgWin)	# Doesn't work with Rgui.exe, but tkwm.deiconify does
}

## Function for the RGui menu
aboutZI <- function (graphical = FALSE)
{
	msg <- getTemp("ZIverstring")
	### TODO: add more information here (copyright, authors, ...)
	if (isTRUE(as.logical(graphical))) {
		dlgMessage(message = msg, title = "About...", icon = "info",
			type = "ok")
	} else cat(msg, "\n")
}

exitZI <- function ()
{
	## This is useful to allow updating the package!
	detach("package:zooimage", unload = TRUE)
	message("zooimage package unloaded; To restart it, issue:\n> library(zooimage)")
}

## Functions for the assistant menu
closeAssistant <- function ()
{
	try(menuDel(getTemp("ZIname")), silent = TRUE)
	try(menuDel("Analyze"), silent = TRUE)
	try(menuDel("Utilities"), silent = TRUE)
	## Destroy the ZooImage Tk window, if it is currently displayed
	#tkWinDel("ZIDlgWin")
	## Eliminate the function to update the progressmeter in that window
	#assignTemp(".progress", list())
}

closeZooImage <- function ()
{
	closeAssistant()
	exitZI()
}

viewManual <- function ()
{
	manual <- file.path(getTemp("ZIetc"), "ZooImageManual.pdf")
	pdfviewer <- getOption( "pdfviewer" )
	if (!is.null(pdfviewer)) {
		if (.Platform$OS.type == "windows") {
            shell.exec(manual)
        } else {
			system(paste(shQuote(getOption("pdfviewer")), shQuote(manual)),
				wait = FALSE)
		}
	} else browseURL(manual)
}

focusR <- function ()
{
	## Switch the focus to the R console
	### TODO: notify this command is not available elsewhere (inactivate menu?)
	if (isRgui()) bringToTop(-1) else
		stop("Not implemented in this environment")
}

focusGraph <- function ()
{
	## Focus to the active R graph (create one if there is no graph device)
	### TODO: notify this command is not available elsewhere (inactivate menu?)
	if (is.null(dev.list())) {
		device <- match.fun(getOption("device"))
		device()
	} else {
		## Activate current graph window
		if (isRgui()) bringToTop() else
			stop("Not implemented in this environment")
	}
}

## Show an assitant dialog box allowing to choose between VueScan and a different
## acquisition program... remember that setting in the registry under Windows
acquireImg <- function ()
{
	## First read the registry to determine which software in recorded there...
 	Asoft <- getOption("ZI.AcquisitionSoftware", "VueScan")
	if (Asoft == "VueScan") {
		opts <- c("VueScan", "Another software...")
		othersoft <- ""
		defval <- "VueScan"
	} else {
		othersoft <- Asoft
       	defval <- basename(othersoft)
		opts <- c("VueScan", defval, "Another software...")
	}

	## Then, show the dialog box
 	#res <- modalAssistant(paste(getTemp("ZIname"), "picture acquisition"),
	#	c("To acquire digital plankton images,",
	#	"you can use a specialized equipment or",
	#	"a digital camera on top of a binocular or",
	#	"a flatbed scanner, ...",
	#	"",
	#	"To pilot a scanner or rework RAW digicam images",
	#	"you can use 'Vuescan'.",
	#	"You can also specify to use your own software.",
	#	"", "", "Use:", ""), init = defval,
	#	options = opts, help.topic = "acquireImg")
	## Analyze result
	#if (res == "ID_CANCEL") return(invisible())
	res <- dlgList(opts, preselect = defval, multiple = FALSE,
		title = "Acquire images with:")$res	
	if (!length(res)) return(invisible())	
	## Did we selected "Another software..."?
	if (res == "Another software...") {
		## Ask for selecting this software
		Asoft <- dlgOpen(title = "Select a program...", multiple = FALSE)$res
		if (!length(Asoft)) return(invisible(NULL)) # Cancelled dialog box
	}
	## Did we selected "VueScan"
	if (res == "VueScan") {
		startPgm("VueScan", switchdir = TRUE)
		options(ZI.AcquisitionSoftware = "VueScan")
		return(invisible(NULL))
	}
	## We should have selected a custom software...
	if (!file.exists(Asoft))
		stop("Program '", Asoft, "' not found!")
	## Start the program
	system(paste('"', Asoft, '"', sep = ""), wait = FALSE)
	## Record it in the registry key
    options(ZI.AcquisitionSoftware = Asoft)
}

importImg <- function ()
{
	# Import images... basically, you can select a series of images in a
	# directory, and the program asks for writing the associated .zim files,
	# or you can access other processes that automatically build .zim files
	# and/or import images/data, including custom processes defined in
	# separate 'ZIEimport' objects (see FlowCAM import routine for an example)
	# Get a list of 'ZIEimport' objects currently loaded in memory

	Images <- selectFile("Img", multiple = TRUE, quote = FALSE,
		title = "Select data to import...")

	## Look if there is at least one image selected
	if (!length(Images)) return(invisible(FALSE))
	dir <- dirname(Images[1])
	Images <- basename(Images)

	has <- function (file, pattern)
		grepl(pattern, Images[1])

	## Determine which kind of data it is
	if (has(Images[1], pattern = "^Import_.*[.]zie$")) {
		if (length(Images) > 1)
			warning("you cannot select more than one ZIE file; using the first one")
		
		return(invisible(zieMake(path = dir, Filemap = Images[1], check = TRUE)))
    
	} else if (has(Images[1], "[.]txt$")) {
		## Special Case for FlowCAM images
		if (length(Images) > 1)
			warning("you cannot select more than one TXT file; using the first one")
		
		## I also need the "ImportTemplate.zie" file in the same path
		txtFile <- Images
		zieTemplate <- file.path(dirname(txtFile), "ImportTemplate.zie")
		if (!checkFileExists(zieTemplate, "zie", force.file = TRUE))
			return(invisible(FALSE))
		
		## Create .zim files + FitVisParameters.csv file for the FlowCAM
		message("Creating .zim files and FitVisParameters.csv...")
		res <- zieCompileFlowCAM(path = dirname(txtFile), Tablefile = txtFile,
			Template = zieTemplate, check.names = FALSE)
		return(invisible(res))
	
	} else if (has(".tif")) {
		pattern <- extensionPattern(".tif")
	
	} else if (has("jpg")) {
        pattern <- extensionPattern("jpg")
	
	} else {
		warning("unrecognized data type!")
		return(invisible(FALSE))
	}

	## If there is no special treatment, just make all required .zim files
	## for currently selected images
	invisible(zimMake(dir = dir, pattern = pattern, images = Images))
}

## TODO: the text appears only on one line on the Mac???
processImg <- function ()
{
	## Display a dialog box telling how to process images using ImageJ
	## When the user clicks on 'OK', ImageJ is started... + the checkbox 'close R'
	#res <- modalAssistant(paste(getTemp("ZIname"), "picture processing"),
	#	c(paste("Once images are acquired and imported into", getTemp("ZIname")),
	#	"(they have correct associated metadata), they must be",
	#	"processed.",
	#	"",
	#	"To do so, start 'ImageJ' (just click 'OK') and select",
	#	paste("the method for your images in 'Plugins -> ", getTemp("ZIname"),
	#		"'.", sep = ""),
	#	"",
	#	"For very large images, or on computers with limited",
	#	"RAM memory, it is advised to close all other programs.",
	#	"Check the option below to close R in this case.", "", ""),
	#	init = "0", check = "Close R before running 'ImageJ'",
	#	help.topic = "processImg")
	## Analyze result
	#if (res == "ID_CANCEL") return(invisible())
	res <- dlgMessage(paste("You will switch now to ImageJ to process your",
		"images. Do you want to continue?"), type = "okcancel")$res
	if (res == "cancel") return(invisible(NULL))
 	## Start ImageJ
	if (!is.null(getOption("ImageEditor")))
		startPgm("ImageEditor", switchdir = TRUE, iconize = TRUE)
	## Do we have to close R?
	#if (res == "1") q()
}

makeZid <- function ()
{
	## Create ZID files, possibly processing imqges first
	## TODO: get the list of all available processes
	## and select it automatically from the ZIM file
	defval <- "Scanner_Gray16"
	
	## Calls the class org.sciviews.zooimage.ZooImageProcessList to get 
	## the list of available processes
	getProcessList <- function () {
		cmd <- sprintf('java -cp .:"%s":"%s" org.sciviews.zooimage.ZooImageProcessList', 
			system.file("imagej", "ij.jar", package = "zooimage"),
			system.file("imagej", "plugins", "_zooimage.jar",
			package = "zooimage"))
		system(cmd , intern = TRUE)
	}	
	processes <- getProcessList()
	opts <- c( processes, "-- None --")
	## Then, show the dialog box
 	#plugin <- modalAssistant(paste(getTemp("ZIname"), "process images"),
	#	c("Process images with associated metadata (ZIM files)",
	#	"in batch mode from one directory and make ZID files.",
	#	"", "Select an image processor:", ""), init = defval,
	#	options = opts, help.topic = "processIJ")
	## Analyze result
	#if (plugin == "ID_CANCEL") return(invisible())
	plugin <- dlgList(opts, preselect = defval, multiple = FALSE,
		title = "Select a batch image processor:")$res	
	if (!length(plugin)) return(invisible(NULL))
	## Select zim file or directory
	dir <- dlgDir()$res
	if (!length(dir)) return(invisible(NULL))
	## Do we need to process the images with ImageJ?
	if (plugin != "-- None --") {
		ijplugin <- function (zimfile, ij.plugin = c("Scanner_Gray16",
			"MacroPhoto_Gray16", "Scanner_Color", "Microscope_Color")){
			ij.plugin <- match.arg(ij.plugin)
			cmd <- sprintf('java -Xmx900m -cp .:"%s":"%s" org.sciviews.zooimage.ZooImage %s "%s"',
				system.file("imagej", "ij.jar", package = "zooimage"),
				system.file("imagej", "plugins", "_zooimage.jar",
					package = "zooimage"), ij.plugin,
					tools::file_path_as_absolute(zimfile))
			return(invisible(system(cmd, intern = TRUE)))
		}		
		## TODO: update a progress bar from ImageJ (using sockets ?)
		ijplugin(dir, ij.plugin = plugin) 
	}

	## Finalize .zid files (and possibly also .zip files by updating their comment)
#    res <- modalAssistant(paste(getTemp("ZIname"), "data processing"),
#		c("You should have processed all your images now.",
#		"The next step is to finalize the .zid files (ZooImage",
#		"Data files). There will be one data file per sample and",
#		"it is all you need for the next part of your work...",
#		"",
#		"Once this step succeed, you can free disk space by",
#		"transferring all files from the _raw subdirectory to",
#		"archives, for instance, DVDs (Apps -> CD-DVD burner).",
#		"",
#        "Warning: the whole _work subdirectory with intermediary",
#		"images will be deleted, and all .zim files will be",
#		"moved to the _raw subdirectory.",
#		"At the end, you should have only .zid files remaining",
#		"in your working directory.", "",
#		"Click 'OK' to proceed (select working directory)...", ""),
#		init = "1", check = "Check vignettes", help.topic = "makeZid")
#	# Analyze result
#	if (res == "ID_CANCEL") return(invisible())
#	# Confirm the directory to process...
#	dir <- dlgDir()$res
#	if (length(dir) == 0) return(invisible())
 	## Do we check the vignettes (only if images were not processed)?
 	check.vignettes <- (plugin == "-- None --")
	## Make .zid files
    cat("\n")
	## TODO: combine the log from ImageJ with this one!
	zidCompressAll(path = dir, check.vignettes = check.vignettes,
		replace = TRUE, delete.source = TRUE)
}

makeTrain <- function ()
{
	## Select samples, and a grouping template... and prepare
	## for making a training set
    ## First read the registry to determine which grouping in recorded there...
 	Grp <- getOption("ZI.DefaultGrouping", "[Basic]")
	## Does this point to an actual file?
	if (file.exists(Grp)) {
		defval <- basename(Grp)
		opts <- c("Basic", "Detailed", "Very_detailed", defval, "Another config...")
		otherGrp <- Grp
	} else {
		defval <- sub("^[[](.+)[]]$", "\\1", Grp)
		opts <- c("Basic", "Detailed", "Very_detailed", "Another config...")
		otherGrp <- ""
	}
	## Then, show the dialog box
 	#res <- modalAssistant(paste(getTemp("ZIname"), "prepare training set"),
	#	c("This step prepares a directory in the hard disk",
	#	"where you will have the opportunity to manually",
	#	"classify vignettes in as many taxa as you want.",
	#	"The hierarchy of the folders and subfolders can",
	#	"be used to represent various levels of classification",
	#	"that the software will be able to use subsequently.",
	#	"",
	#	"You must specify: (1) a grouping scheme to start with,",
	#	"(2) a base directory where to locate the training set,",
	#	"(3) a series of .zid files as source of vignettes.", "",
	#	"Use the following grouping scheme:", ""), init = defval,
	#	options = opts, help.topic = "makeTrain")

	## Analyze result
	#if (res == "ID_CANCEL") return(invisible())
	res <- dlgList(opts, preselect = defval, multiple = FALSE,
		title = "Select the default classes to use to initialize your training set:")$res	
	if (!length(res)) return(invisible(NULL))

	## Did we selected "Another config..."?
	if (res == "Another config...") {
		## Ask for selecting a .zic file containing the config
        otherGrp <- selectFile("Zic", multiple = FALSE, quote = FALSE,
			title = "Select a .zic file...")
		if (!length(otherGrp)) return(invisible(NULL))
		## Cancelled dialog box
		res <- otherGrp
	} else if (res %in% c("Basic", "Detailed", "Very_detailed")) {
		## Did we selected a standard scheme?
		res <- paste("[", res, "]", sep = "")
	} else res <- Grp  # We should have selected the previously recorded scheme...

	## Save this config for later use
    options(ZI.DefaultGrouping = res)

	## Ask for the base directory
    dir <- dlgDir()$res
	if (!length(dir)) return(invisible(NULL))

	## Ask for a subdir for this training set
	subdir <- dlgInput("Subdirectory where to create the training set:",
		default = "_train")$res
	if (!length(subdir)) return(invisible(NULL))

	## Ask for the .zid files
    zidfiles <- selectFile(type = "Zid", multiple = TRUE, quote = FALSE)
	if (!length(zidfiles)) return(invisible(NULL))

	## Prepare the training set
	prepareTrain(file.path(dir, subdir), zidfiles, template = res)
	imageViewer(file.path(dir, subdir, "_"))

	## Remember the directory...
	assignTemp("ZI.TrainDir", file.path(dir, subdir))
}

## Read a training set and create a ZITrain object
collectTrain <- function ()
{
	## Get a possibly saved directory as default one
	dir <- getTemp("ZI.TrainDir")
	if (is.null(dir) || !file.exists(dir) || !file.info(dir)$isdir)
		dir <- getwd()
	## Ask for a base directory of a training set...
	dir <- dlgDir(default = dir, title = paste("Select a", getTemp("ZIname"),
		"training set base dir"))$res
	if (!length(dir) || !file.exists(dir) || !file.info(dir)$isdir)
		return(invisible(NULL))
	
	## Ask for a name for this ZITrain object
	name <- dlgInput("Name for the ZITrain object to create in the global environment:",
		default = "ZItrain")$res
	if (!length(name)) return(invisible(FALSE))
	name <- make.names(name)	# Make sure it is a valid name!
	
	## Get the training set and save it in .GlobalEnv under the provided name
	res <- getTrain(dir, creator = NULL, desc = NULL, keep_ = FALSE)
	.assignGlobal(name, res)
	
	## Remember the object name
	assignTemp("ZI.TrainName", name)
	
	## Print informations about this training set
	message("Manual training set data collected in '", name, "'")
	cat("\nClassification stats:\n")
	print(table(res$Class))
	cat("\nProportions per class:\n")
	print(table(res$Class) / length(res$Class) * 100)
}

## Add data to an existing training set
addVigsToTrain <- function ()
{
	## Select zid or zidb files to add in the training set
	zidb <- selectFile(type = "ZidZidb", multiple = TRUE, quote = FALSE)
	if (!length(zidb)) return(invisible(NULL))
	
	## Select the training set in which we add new vignettes
	dir <- getTemp("ZI.TrainDir")
	if (is.null(dir) || !file.exists(dir) || !file.info(dir)$isdir)
		dir <- getwd()
	## Ask for a base directory of a training set...
	dir <- dlgDir(default = dir, title = paste("Select a", getTemp("ZIname"),
		"training set base dir"))$res
	if (!length(dir) || !file.exists(dir) || !file.info(dir)$isdir)
		return(invisible(NULL))
	
	## Extract vignettes in the training set in a _NewVignettesX directory
	message("Adding vignettes from these files to _ subdir...")
	addToTrain(traindir = dir, zidbfiles = zidb)
}

## New version to accept variables selection and/or new formula 1.2-2
## TODO: avoid duplication of code here
makeClass <- function ()
{
 	## Create a classifier, using a ZI1Class object (new version)
	## Ask for an algorithm + additional parameters
	## Return a ZIClass object
	defval <- "linear discriminant analysis"
	opts <- c("linear discriminant analysis",
			  "recursive partitioning (tree)",
			  "k-nearest neighbour",
			  "learning vector quantization",
			  "neural network",
			  "random forest",
			  "Variables Selection")

 	#res <- modalAssistant(paste(getTemp("ZIname"), "make classifier"),
	#	c("This is a simplified version of the classifiers",
	#	"where you just need to select one algorithm.",
	#	"Warning! Many algorithms have parameters to be",
	#	"fine-tuned before efficient use... and this must be",
	#	"done for each specific data set! Here, only default",
	#	"parameters that have proven efficient with plankton",
	#	"are applied automatically. Some methods already work",
	#	"pretty well that way.",
	#	"", "Learn using an algorithm:", ""), init = defval,
	#	options = opts, help.topic = "makeClass")
	## Analyze result
	#if (res == "ID_CANCEL") return(invisible())
	res <- dlgList(opts, preselect = defval, multiple = FALSE,
		title = "Select an algorithm for creating your classifier:")$res	
	if (!length(res)) return(invisible(NULL))

	if (res != "Variables Selection") {
		## Use default values for the classifier creation
		warnings("These defaults variables are used : logArea, Mean, StdDev, ",
			"Mode, Min, Max, logPerim., logMajor, logMinor, Circ., logFeret, ",
			"IntDen, Elongation, CentBoxD, GrayCentBoxD, CentroidsD, Range, ",
			"MeanPos, SDNorm, CV")
		## Compute algorithm from res
		algorithm <- switch(res,
			`linear discriminant analysis` = "lda",
			`recursive partitioning (tree)` = "rpart",
			`random forest` = "randomForest",
			`support vector machine` = "svm",
			`k-nearest neighbour` = "ipredknn",
			`learning vector quantization` = "mlLvq",
			`neural network` = "mlNnet")

		## Look if we have a manual training set object defined
		ZIT <- getTemp("ZI.TrainName")
		if (is.null(ZIT)) ZIT <- ""
		## Ask for a ZITrain object
		ZIT <- selectObject("ZITrain", multiple = FALSE, default = ZIT,
			title = "Choose one ZITrain objects:")
		if (!length(ZIT) || (length(ZIT) == 1 && ZIT == ""))
			return(invisible(NULL))
		## Ask for a name for this ZIClass object
		name <- dlgInput("Name for the ZIClass object to create in the global environment:",
			default = "ZIclass")$res
		if (!length(name)) return(invisible(NULL))
		name <- make.names(name)	# Make sure it is a valid name!
		## Calculate results
		res <- ZIClass(data = get(ZIT, envir = .GlobalEnv), algorithm = algorithm)
	} else {
		## Options if 'Variables Selection is selected v 1.2-2
		opts <- c("linear discriminant analysis",
				"recursive partitioning (tree)",
				"k-nearest neighbour",
				"learning vector quantization",
				"neural network",
				"random forest")
		## Dialog box if 'Variables Selection' is selected v1.2-2
		#res <- modalAssistant(paste(getTemp("ZIname"), "make classifier"),
		#	c("This is a simplified version of the classifiers",
		#	"where you just need to select one algorithm.",
		#	"Warning! Many algorithms have parameters to be",
		#	"fine-tuned before efficient use... and this must be",
		#	"done for each specific data set!",
		#	"",
		#	"Here, you can select",
		#	"variables to use for the classifier creation.",
		#	"",
		#	"Warning! Select only pertinent and useful measurements.",
		#	"", "Learn using an algorithm:", ""), init = defval,
		#	options = opts, help.topic = "makeClass")
		#if (res == "ID_CANCEL") return(invisible())
		res <- dlgList(opts, preselect = defval, multiple = FALSE,
			title = "Select an algorithm for creating your classifier:")$res	
		if (!length(res)) return(invisible(NULL))

		## Compute algorithm from res
		algorithm <- switch(res,
			`linear discriminant analysis` = "lda",
			`recursive partitioning (tree)` = "rpart",
			`random forest` = "randomForest",
			`support vector machine` = "svm",
			`k-nearest neighbour` = "ipredknn",
			`learning vector quantization` = "mlLvq",
			`neural network` = "mlNnet")

		## Look if we have a manual training set object defined
		ZIT <- getTemp("ZI.TrainName")
		if (is.null(ZIT)) ZIT <- ""
		## Ask for a ZITrain object
		ZIT <- selectObject("ZITrain", multiple = FALSE, default = ZIT,
			title = "Choose one ZITrain objects:")
		if (length(ZIT) == 0 || (length(ZIT) == 1 && ZIT == ""))
			return(invisible(NULL))
		## Ask for a name for this ZIClass object
		name <- dlgInput("Name for the ZIClass object to create:",
			title = "Creating a classifier", default = "ZIclass")$res
		if (!length(name)) return(invisible(NULL))
		name <- make.names(name)	# Make sure it is a valid name!
		## Calculate formula using variables of the training set

### TODO: change this: do not return a formula, but a list of variables to select + "Class"!
		## variables selection for the classifier
		## TODO: select variables to use for classification instead of returning a formula!
		selectVars <- function (ZITrain,
		calc.vars = getOption("ZI.calcVars", "calcVars")) {
			## ZITrain must be a ZItrain object
			if (!inherits(ZITrain, "ZITrain")) {
				warning("'ZITrain' must be a 'ZITrain' object")
				return(character(0))
			}
		
			calcfun <- match.fun(as.character(calc.vars)[1])
			## Parameters measured on particles and new variables calculated
			mes <- as.vector(colnames(calcfun(ZITrain)))
			presel <- c("Id", "FIT_Cal_Const", "Item", "FIT_Raw_Area",
				"FIT_Raw_Feret_Max", "FIT_Raw_Feret_Min", "FIT_Raw_Feret_Mean",
				"FIT_Raw_Perim", "FIT_Raw_Convex_Perim", "FIT_Feret_Max_Angle",
				"FIT_Feret_Min_Angle", "FIT_Avg_Red", "FIT_Avg_Green", "FIT_Avg_Blue",
				"FIT_PPC", "FIT_Ch3_Peak", "FIT_Ch3_TOF", "FIT_Ch4_Peak", "FIT_Ch4_TOF",
				"FIT_SaveX", "FIT_SaveY", "FIT_PixelW", "FIT_PixelH", "FIT_CaptureX",
				"FIT_CaptureY", "FIT_Edge_Gradient", "FIT_Timestamp1", "FIT_Timestamp2",
				"FIT_Source_Image", "FIT_Calibration_Image", "FIT_High_U32",
				"FIT_Low_U32", "FIT_Total", "FIT_Red_Green_Ratio",
				"FIT_Blue_Green_Ratio", "FIT_Red_Blue_Ratio", "FIT_Ch2_Ch1_Ratio",
				"X.Item.1", "X", "Y", "XM", "YM", "BX", "BY", "Width", "Height",
				"Angle", "XStart", "YStart", "Count",  "Label", "Dil", "Class")
			DontKeep <-  dlgList(mes, preselect = presel, multiple = TRUE,
				title = "Select variable you don't want to use in the classification")$res
			
			## Selection of features for the creation of the classifier
		#	keep <- dlgList(mes, preselect = c("ECD", "FIT_Area_ABD",
		#		"FIT_Diameter_ABD", "FIT_Volume_ABD", "FIT_Diameter_ESD",
		#		"FIT_Volume_ESD", "FIT_Length", "FIT_Width", "FIT_Aspect_Ratio",
		#		"FIT_Transparency", "FIT_Intensity", "FIT_Sigma_Intensity",
		#		"FIT_Sum_Intensity", "FIT_Compactness", "FIT_Elongation",
		#		"FIT_Perimeter", "FIT_Convex_Perimeter", "FIT_Roughness",
		#		"FIT_Ch1_Peak", "FIT_Ch1_TOF", "FIT_Ch2_Peak", "FIT_Ch2_TOF",
		#		"Area", "Mean", "StdDev", "Mode", "Min", "Max", "Perim.", "Width",
		#		"Height", "Major", "Minor", "Circ.", "Feret", "IntDen", "Median",
		#		"Skew", "Kurt", "Elongation", "CentBoxD", "GrayCentBoxD", "CentroidsD",
		#		"Range", "MeanPos", "SDNorm", "CV", "logArea", "logPerim.", "logMajor",
		#		"logMinor", "logFeret"),
		#		multiple = TRUE, title = "Select variables to use for classification")$res
			
			## Creation of one formula for classifier calculation
			keep <- mes[!mes %in% DontKeep]
			res <- as.formula(paste("Class ~ ", paste(keep, collapse = "+")))
			return(res)
		}
		form <- selectVars(get(ZIT, envir = .GlobalEnv, inherits = FALSE))
		## Calculate results using formula created by variables selection
		res <- ZIClass(form, data = get(ZIT, envir = .GlobalEnv), algorithm = algorithm)
	}
	## Store the resulting object
	.assignGlobal(name, res)
	## Print results
	print(res)
	cat("\n")
	## Remember that ZIClass object
    assignTemp("ZI.ClassName", name)
}

## Analyze confusion matrix
analyzeClass <- function ()
{
	## Analyze a classifier, using a ZI1Class object (new version)
	## Ask for an option of analysis
 	defval <- "Print Confusion Matrix"
	opts <- c("Print Confusion Matrix", "Summarize", "Plot Confusion Matrix",
		"Plot F-score", "Plot Dendrogram", "Plot Precision/recall")
	## Then, show the dialog box
 	#res <- modalAssistant(paste(getTemp("ZIClass"), "Analyze a classifier"),
	#	c("This is a simplified version of the analysis of classifiers",
	#	"where you just need to select one classifier.",
	#	"These options provide some tools to analyze your classifers.",
	#	"", "Select a classifer and a tool:", ""), init = defval,
	#	options = opts)
	## Analyze result
	#if (res == "ID_CANCEL") return(invisible()) # not error message is 'cancel'
	res <- dlgList(opts, preselect = defval, multiple = FALSE,
		title = "Select a classifier to be analyzed:")$res	
	if (!length(res)) return(invisible(NULL))
		
 	## Analyze a classifier... currently, only calculate the confusion matrix
	## and edit it
	ZIC <- selectObject("ZIClass", multiple = FALSE,
		title = "Choose one ZIClass object:")
	if (!length(ZIC))
		stop("No classifier. Please, create one first!")

	ZIC <- get(ZIC, envir = .GlobalEnv)
	conf <- confusion(ZIC)
	switch(res,
		`Print Confusion Matrix` = print(conf),
		`Summarize` = print(summary(conf)),
		`Plot Confusion Matrix` = plot(conf, type = "image"),
		`Plot F-score` = plot(conf, type = "barplot"),
		`Plot Dendrogram` = plot(conf, type = "dendrogram"),
		`Plot Precision/recall` = plot(conf, type = "stars"))
	return(invisible(conf))
}

## Extract vignettes from zid files to respective directories
vignettesClass <- function ()
{
	## Ask for the base directory
    defdir <- getTemp("ZI.BaseDir", default = getwd())
	basedir <- dlgDir(default = defdir,
		title = "Select the base directory for the test set")$res
	if (!length(basedir)) return(invisible(NULL))

	## Ask for a subdir for this training set
	subdir <- dlgInput("Subdirectory where to create the test set:",
		default = "_test")$res
	if (!length(subdir)) return(invisible(NULL))
	testdir <- file.path(basedir, subdir)
	if (file.exists(testdir))
		stop("The directory '", testdir,
			"' already exists! Please, restart and specify a new one")
	
	## Select .zid files to be classified
	zid <- selectFile(type = "ZidZidb", multiple = TRUE, quote = FALSE)
	if (!length(zid)) return(invisible(NULL))
	
	## Look if we have a classifier object defined
	zic <- getTemp("ZI.ClassName", default = "")
	zic <- selectObject("ZIClass", multiple = FALSE, default = zic,
		title = "Choose a classifier (ZIClass object):")
	if (!length(zic)) return(invisible(FALSE))
	## Save this choice for later reuse
	assignTemp("ZI.ClassName", zic, replace.existing = TRUE)
	zicObj <- get(zic, envir = .GlobalEnv)

	## Sort vignettes in the different directories, as predicted by the classifier
	prepareTest(testdir, zid, template = zicObj, classes = zicObj)
	
	## Remember the directory...
	assignTemp("ZI.BaseDir", basedir)
	assignTemp("ZI.TestDir", testdir)
	
	## Explain what to do next...
	message("Vignettes classified in '", testdir, "'")
	message("View them in your favorite file browser (and possibly correct classification manually)")
	
	## Classify vignettes  
#	if (length(zid) > 1) {
#		classVignettesAll(zidfiles = zid, Dir = "_manuValidation",
#			ZIClass = zicObj)
#	} else { # Possibly apply a filter		
#		## Give a name for the final directory
#		finalDir <- dlgInput("Name for the automatic classification directory:",
#			default = noExtension(zid), title = "Parameter filter")$res
#		if (!length(finalDir)) return(invisible(NULL))
#		
#		## Read the zid file
#		ZIDat <- zidDatRead(zid)
#    
#		## Select a parameter to use for the threshold
#		threshold <- createThreshold(ZIDat = ZIDat)     
#		if (length(threshold)) {
#			classVignettes(zidfile = zid, Dir = finalDir,ZIClass = zicObj,
#				ZIDat = ZIDat, Filter = threshold)
#		} else {
#			classVignettes(zidfile = zid, Dir = finalDir, ZIClass = zicObj)
#		}
#	}
}

## Edit a samples description file... or create a new one!
editDescription <- function ()
{
	#res <- modalAssistant(paste(getTemp("ZIname"), "edit samples description"),
	#	c("Samples are about to be analyzed and collected together",
	#	"to form a series.",
	#	paste(getTemp("ZIname"), "needs to know which samples should be"),
	#	"collected into the same series and you must provide",
	#	"metadata information (especially date and time of",
	#	"collection, location of the sampling stations, or",
	#	"possibly temperature, salinity, turbidity, etc. that",
	#	"were recorded at the same time as these samples).",
	#	"",
	#	"A .zis file (by default, Description.zis) needs to be",
	#	"created and edited for each of the considered series.",
	#	"You can here edit, or create a new samples description",
	#	"file from the template.", "",
	#	"Click 'OK' to edit a samples description file now...", ""),
	#	init = "1", check = "New description file from template.",
	#	help.topic = "editDescription")
	## Analyze result
	#if (res == "ID_CANCEL") return(invisible())
	res <- dlgMessage(paste("Create a new description file from scratch?"),
		type = "yesnocancel")$res
	if (res == "cancel") return(invisible(NULL))
	## Edit/create the description file...
	if (res == "yes") {	# Create a Zis file ()take care: was "1" for modalAssistant!
		res <- dlgSave(default = "Description.zis",
			title = "Create a new ZIS file",
			filters = matrix(c("ZooImage samples description", ".zis"),
			ncol = 2, byrow = TRUE))$res
		if (!length(res)) return(invisible(NULL))
		if (regexpr("[.][zZ][iI][sS]$", res) < 0) res <- paste(res, ".zis",
			sep = "")
		zisfile <- zisCreate(res)
	} else { # Edit a Zis file
	    zisfile <- zisEdit(NULL)
	}
	## Remember the last zis file
    assignTemp("ZI.LastZIS", zisfile)
}

processSamples <- function()
{
	## Ask for a description.zis file, look at all samples described there
	## Calculate abundances, total and partial size spectra and possibly biomasses
	## Get the last edited description.zis file
	## Get a possibly saved directory as default one
	zisfile <- getTemp("ZI.LastZIS")
	if (is.null(zisfile) || !file.exists(zisfile))
		zisfile <- ""
	## Ask for a file
	if (zisfile != "") {	
		zisfile <- dlgOpen(default = zisfile, title = "Select a ZIS file",
			filters = matrix(c("ZooImage samples description", ".zis"),
			ncol = 2, byrow = TRUE))$res	
	} else if (file.exists(file.path(getwd(), "Description.zis"))) {
		zisfile <- dlgOpen(default = file.path(getwd(), "Description.zis"),
			title = "Select a ZIS file",
			filters = matrix(c("ZooImage samples description", ".zis"),
			ncol = 2, byrow = TRUE))$res	
	} else {
		zisfile <- dlgOpen(title = "Select a ZIS file",
			filters = matrix(c("ZooImage samples description", ".zis"),
			ncol = 2, byrow = TRUE))$res	
	}
	if (!length(zisfile)) return(invisible(NULL))

	## Add Kevin to use manual validation 2010-08-03
	## Option dialog box
	#res <- modalAssistant(paste(getTemp("ZIname"), "samples processing"),
	#	c(
	#		"Each sample registered in the description.zis file",
	#		"will be processed in turn to extract ecological",
	#		"parameters (abundances, biomasses, size spectra).",
	#		"",
	#		"If you want to save calculation done on each",
	#		"particle individually, check the option below.",
	#		"",
	#		"Click 'OK' to proceed...", ""
	#	), init = "0",
	#	options = "Manual Validation", check = "Save individual calculations",
	#	help.topic = "processSamples")
	## Analyze result
	#if (res == "ID_CANCEL") return(invisible())
	res <- dlgMessage(paste("Save also calculations done on each particle individually?"),
		type = "yesnocancel")$res
	if (res == "cancel") return(invisible(NULL))
	## Do we save individual calculations?
	if (res == "yes")	# Note that for modalAssistant, it was "1"!
		exportdir <- dirname(zisfile) else exportdir <- NULL
	## Added by Kevin for semi auto classif
	## Do we use Semi automatic classification?
	if (res == "Manual Validation") {
		#res <- modalAssistant(paste(getTemp("ZIname"), "samples processing"),
		#c(
		#	"Each sample registered in the description.zis file",
		#	"will be processed in turn to extract ecological",
		#	"parameters (abundances, biomasses, size spectra)",
		#	"after manual validation of automatic predictions",
		#	"done in the '_manualValidation' directory", 
		#	"",
		#	"If you want to save calculation done on each",
		#	"particle individually, check the option below.",
		#	"",
		#	"Click 'OK' to proceed...", ""
		#), init = "0",
		#check = "Save individual calculations", help.topic = "processSamples")
		## Analyze result
		#if (res == "ID_CANCEL") return(invisible())
		res <- dlgMessage(paste("Save also calculations done on each particle individually?"),
			type = "yesnocancel")$res
		if (res == "cancel") return(invisible(NULL))
		## Do we save individual calculations?
		if (res == "yes") # Note that for modalAsisstant, it was "1"!
			exportdir <- dirname(zisfile) else exportdir <- NULL
		
## TODO: change this!
#		## Select the directory where manual validation is done
#		dir <- getTemp("ZI.TrainDir")
#		if (is.null(dir) || !file.exists(dir) || !file.info(dir)$isdir)
#			dir <- getwd()
#		## Ask for a base directory of a training set...
#		dir <- dlgDir(default = dir, title = paste("Select a",
#			getTemp("ZIname"), "Manual validation base dir"))$res
#		if (!length(dir) || !file.exists(dir) || !file.info(dir)$isdir)
#			return(invisible(NULL))
#		## Read the directory
#		ZIManTable <- ZIManRead(dir)
#		message("Read the manual validation directory...\n-- Done --")		
#		ManValid <- TRUE
#	} else {
#		## Classification without any manual validation
		ManValid <- FALSE
	} 
	
	## Get a list of samples from the description file
	smpdesc <- zisRead(zisfile)
	smplist <- listSamples(smpdesc)
	if (!length(smplist) || smplist == "")
		stop("No sample found in the description file!")
	
	## Are there corresponding .zid files for all samples?
	zisdir <- dirname(zisfile)
	if (zisdir == ".") zisdir <- getwd()
	zidfiles <- file.path(zisdir, paste(smplist, ".zid", sep = ""))
	if (!all(file.exists(zidfiles)) ||
		!all(regexpr("[.][zZ][iI][dD]$", zidfiles) > 0))
		stop("One or more .zid files do not exist or is invalid!")
	
	## Get a classifier
	ZIC <- getTemp("ZI.ClassName")
	if (is.null(ZIC)) ZIC <- ""
	ZIC <- selectObject("ZIClass", multiple = FALSE, default = ZIC,
		title = "Choose a classifier (ZIClass object):")
	if (!length(ZIC) || (length(ZIC) == 1 && ZIC == ""))
		return(invisible(NULL))
	ZICobj <- get(ZIC, envir = .GlobalEnv)	
	
	## Read a conversion table from disk (from /etc/Conversion.txt)
	## or an other position
	## First read the options to determine which file to use...
	ConvFile <- getOption("ZI.ConversionFile", file.path(getTemp("ZIetc"),
		"Conversion.txt"))
	## Does this file exists?
	if (!file.exists(ConvFile) || ConvFile == "")
		ConvFile <- file.path(getTemp("ZIetc"), "Conversion.txt")
	## Ask for selecting a Conversion file
	ConvFile2 <- dlgOpen(default = ConvFile,
		title = "Select a conversion file...", multiple = FALSE,
		filters = matrix(c("Biomass Conversion table (*Conversion.txt)", "Conversion.txt"),
		ncol = 2, byrow = TRUE))$res
	if (!length(ConvFile2)) return(invisible(NULL)) # Cancelled dialog box
	
	## Read the data from this table
	conv <- read.table(ConvFile2, header = TRUE, sep = "\t")
	
	## Save this config for later use
	options(ZI.ConversionFile = ConvFile2)
	
	## Get class breaks for size spectra
	brks <- dlgInput("Breaks for size spectrum classes (empty for no spectrum):",
		default = "seq(0.25, 2, by = 0.1)")$res
 	if (!length(brks)) return(invisible(NULL))
	brks <- eval(parse(text = brks))

	## Get a name for the variable containing results
	name <- dlgInput("Name for the ZIRes object to create in the global environment:",
		default = "ZIres")$res
	if (!length(name)) return(invisible(NULL))
	name <- make.names(name)
	## Add Kevin for manual validation
	if (!isTRUE(as.logical(ManValid))) ZIManTable <- NULL 
	## TODO: we need at least keep and detail
	res <- processSampleAll(path = dirname(zisfile), #ZIClass = ...,
		biomass = conv, breaks = brks)
	## TODO: possibly export result in a file...
	
	
	## Assign this result to the variable
	.assignGlobal(name, res)
	## Remember the name of the variable
	assignTemp("ZI.LastRES", name)
}

viewResults <- function ()
{
 	## Make graphic representations of results...
	ZIR <- getTemp("ZI.LastRES")
	if (is.null(ZIR)) ZIR <- ""
	ZIR <- selectObject("ZIRes", multiple = FALSE, default = ZIR,
		title = "Choose one ZIRes object:")
	if (!length(ZIR) || (length(ZIR) == 1 && ZIR == ""))
		return(invisible(NULL))
	## Get the object
	ZIR <- get(ZIR, envir = .GlobalEnv)
	## Ask for selecting items in the list and make these graphs
	## Compute the list of graphs
	vars <- names(ZIR)
	## Eliminate variables that cannot be plotted...
	vars <- vars[-(1:25)]
	vars <- vars[!vars == "Note"]
	## Add the spectra graphs
	spec <- attr(ZIR, "spectrum")
	varspec <- paste("spectrum of", names(spec))
	Vars <- c(vars, varspec)
	Keep <- dlgList(Vars, multiple = TRUE, title = "Select 1 to 12 graphs:")$res
	lKeep <- length(Keep)
	if (lKeep == 0) return(invisible())
	if (lKeep > 12) {
		Keep <- Keep[1:12]
		lKeep <- 12
	}
	## If there are spectrums, also ask for partial spectrums
	if (any(regexpr("^spectrum of ", Keep) + 1)) {
		pspec <- names(spec[[1]])
		## Replace total by [none] in this list
		pspec[pspec == "total"] <- "[none]"
		Pspec <- dlgList(pspec, multiple = FALSE,
			title = "Select taxon for partial spectrum:")$res
		if (!length(Pspec)) return(invisible(NULL))
	} else Pspec <- "[none]"
	## Do the graphs
	## Determine number of rows and columns
	nc <- round((lKeep + 1) / 3)
	if (nc > 3) nc <- 3
	if (lKeep == 8) nc <- 2
	nr <- c(1, 2, 3, 2, 3, 3, 3, 4, 3, 4, 4, 4)[lKeep]
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	par(mfrow = c(nr, nc))
	for (i in 1:lKeep) {
    	## Determine if it is a x/y graph, or a size spectrum
		if (regexpr("^spectrum of ", Keep[i]) > 0) { # Size spectrum
			Ser <- sub("^spectrum of ", "", Keep[i])
			plot(spec[[Ser]][["total"]], lwd = 3, col = "gray", type = "h",
				main = Ser, ylab = "Frequency")
			if (Pspec != "[none]"){
				Spec <- spec[[Ser]][[Pspec]]
				Spec[Spec == 0] <- NA
				points(Spec, lwd = 6, col = 2, type = "h")
			}
		} else { # x/y graph
			 ## If there is NA in a variable, the plot generates an error
			 Xdat <- ZIR[, "Date"]
			 Ydat <- ZIR[, Keep[i]]
			 if (all(is.na(Xdat)) || all(is.na(Ydat))) {
			    plot(0:1, 0:1, type = "n", xlab = "", ylab = "", xaxt = "n",
					yaxt = "n", main = Keep[i])
			    text(0.5, 0.5, "No data!", adj = c(0.4, 0.5))
			} else {
			 	plot(Xdat, Ydat, xlab = "Date", ylab = Keep[i], main = Keep[i])
			}
		}
	}
}

exportResults <- function ()
{
 	## Export one or more ZIRes objects to text files...
    res <- selectObject("ZIRes", multiple = TRUE,
		title = "Choose one or more ZIRes objects:")
	if (!length(res) || (length(res) == 1 && res == ""))
		return(invisible(NULL))
	## Select a directory where to place these files
	dir <- dlgDir()$res
	if (!length(dir)) return(invisible(NULL))
	filenames <- file.path(dir, res)
	## Export them there
	for (i in 1:length(res)) {
    	## We write those tables:
		## 1) Results [name].txt
		## 2) Metadata [name]_metadata.txt
		## 3) Size spectra [name]_spect_[sample].txt
		obj <- get(res[i], envir = .GlobalEnv)
		write.table(obj,  file = paste(filenames[i], "_AbdBio.txt", sep = ""),
			quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = getDec(),
			row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
		spc <- attr(obj, "spectrum")
		spcnames <- names(spc)
		if (!is.null(spcnames) && length(spcnames) > 0) {
			for (j in 1:length(spcnames)) {
				## Construct a data frame
				spc1 <- spc[[spcnames[j]]]
				breaks <- attr(spc1, "breaks")
				breaks <- breaks[1:(length(breaks) - 1)]
				spctab <- as.data.frame(spc1)
				spctab <- spctab[ , seq(2, ncol(spctab), by = 2)]
				names(spctab) <- names(spc1)
				spctab <- data.frame(breaks = breaks, spctab)
				write.table(spctab,
					file = paste(filenames[i], "_Spectrum_", spcnames[j],
					".txt", sep = ""), quote = FALSE, sep = "\t", eol = "\n",
					na = "NA", dec = getDec(), row.names = FALSE,
					col.names = TRUE, qmethod = c("escape", "double"))
			}
		}
	}
	message(i, "ZIRes object(s) exported in'", dir, "'")
}

loadObjects <- function ()
{
	file <- selectFile("RData", multiple = FALSE, quote = FALSE,
		title = "Select a RData file...")
	if (!length(file)) return(invisible(NULL)) # Cancelled dialog box
	if (file.exists(file)) load(file, envir = .GlobalEnv)
}

saveObjects <- function ()
{
	Objects <- selectObject(c("ZIDat", "ZIDesc", "ZITrain", "ZIClass", "ZIRes",
		"ZIRecode"), multiple = TRUE,
		title = paste("Choose", getTemp("ZIname"), "object(s):"))
	if (!length(Objects) || (length(Objects) == 1 && Objects == ""))
		return(invisible(FALSE))
	file <- dlgSave(default = paste(getTemp("ZIname"), ".RData", sep = ""),
		title = paste("Save", getTemp("ZIname"), "data under..."),
		multiple = FALSE, filters = matrix(c("R data", ".RData"),
		ncol = 2, byrow = TRUE))$res
	if (!length(file)) return(invisible(NULL))
	if (regexpr("[.][rR][dD][aA][tT][aA]$", file) < 0)
		file <- paste(file, ".RData", sep = "")
	save(list = Objects, file = file, compress = TRUE)
}

listObjects <- function ()
{
    varlist <- objects(pos = 1)
	if (!length(varlist))
		stop("No objects currently loaded in memory!")
	Filter <- NULL
	for (i in 1:length(varlist)) Filter[i] <- inherits(get(varlist[i]),
		c("ZIDat", "ZIDesc", "ZITrain", "ZIClass", "ZIRes", "ZIRecode"))
	varlist <- varlist[Filter]
	if (!length(varlist)) {
		stop("No ", getTemp("ZIname"), " objects currently loaded in memory!")
	} else {
    	print(varlist)
	}
}

removeObjects <- function ()
{
	Objects <- selectObject(c("ZIDat", "ZIDesc", "ZITrain", "ZIClass", "ZIRes",
		"ZIRecode"), multiple = TRUE,
		title = paste(getTemp("ZIname"), "object(s) to remove:"))
	if (!length(Objects) || (length(Objects) == 1 && Objects == ""))
		return(invisible(FALSE))
	rm(list = Objects, envir = .GlobalEnv)
}

calib <- function ()
{
	## Select calibration file (*.tif or *.pgm) and calculate White/Black point
	file <- selectFile("TifPgm", multiple = FALSE, quote = FALSE,
		title = "Select a calibration image...")
	if (!length(file)) return(invisible(NULL)) # Cancelled
	if (file.exists(file)) {
		message("Calibrating gray scale... [", basename(file), "]")
		flush.console()
		res <- calibrate(file)
		message("WhitePoint=", round(res["WhitePoint"]))
		message("BlackPoint=", round(res["BlackPoint"]))
		if (length(attr(res, "msg")) > 0)
			message("\nTake care:")
		message(paste(attr(res, "msg"), collapse = "\n"))
	}
}

optInOutDecimalSep <- function ()
{
	## Define the default numeric decimal separator for input and output
	Dec <- getDec()
	## Possibly ask for another one
	DecList <- c(".", ",")
	DecSel <- dlgList(DecList, preselect = Dec, multiple = FALSE,
		title = "In/Out decimal separator")$res
	## Is the cancel button pressed, or is it still the same decimal separator
	if (!length(DecSel) || DecSel == Dec) return(invisible(Dec))
	## Record it in options
    options(OutDec = DecSel)
    ## Indicate change
    cat("In/Out decimal separator changed to '", DecSel, "'\n", sep = "")
 	return(invisible(DecSel))
}


###### Not in menus yet! ##################
## Subpart of zid file and return a subtable corresponding to the threshold
## TODO: is this really a top-menu function... or is it supposed to be used elsewhere?
#subpartZIDat <- function ()
#{
#    ## Select files to use
#    zidFile <- selectFile(type = "Zid", multiple = FALSE, quote = FALSE)
#	if (!length(zidFile)) return(invisible(NULL))
#
#    ## Read the zid file
#    zid <- zidDatRead(zidFile)
#
#    ## Select a parameter to use for the threshold
#    threshold <- createThreshold(ZIDat = zid)    
#
#    ## Apply the thresold
#    res <- subpartThreshold(ZIDat = zid, Filter = threshold)
#    return(res)
#}
