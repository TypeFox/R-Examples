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
## along with ZooImage.  If not, see <http://www.gnu.org/licenses/>.

## Loading and unloading ZooImage
.onLoad <- function (libname, pkgname)
{
	if (!interactive()) options(ZIAssistant  = FALSE)

	## Use the SciViews style for dialog boxes
	options(guiStyle = "SciViews")

	## Did we redefined the ZooImage config?
	redef <- getOption("ZI.redefine")
	if (is.null(redef)) redef <- FALSE else redef <- TRUE
	options(ZI.redefine = NULL)

	## Create some strings in TempEnv
	ZIversion <- packageDescription("zooimage", fields = "Version")
	assignTemp("ZIversion", ZIversion)

	ZIname <- getTemp("ZIname")
	if (!redef || is.null(ZIname)) ZIname <- "ZooImage"
	assignTemp("ZIname", ZIname)
	assignTemp("ZIverstring", paste(ZIname, "version", ZIversion))

	ZIetc <- getTemp("ZIetc")
	if (!redef || is.null(ZIetc))
		ZIetc <- system.file("etc", package = "zooimage")
	assignTemp("ZIetc", ZIetc)

	ZIgui <- getTemp("ZIgui")
	if (!redef || is.null(ZIgui))
		ZIgui <- system.file("gui", package = "zooimage")
	assignTemp("ZIgui", ZIgui)

	## Windows specific things
	#if (isWin()) {
		#if (interactive()) {
		#	ZIico <- getTemp("ZIico")
		#	if (!redef || is.null(ZIgui))
		#		ZIico <- tk2ico.create(file.path(getTemp("ZIgui"),
		#			"ZooImage.ico"))
		#	assignTemp("ZIico", ZIico)
		#}

		## Make sure there is a key for ZooImage in the registry
		## PhG: what is the purpose of this code?
		#ZIkey <- "HKEY_LOCAL_MACHINE\\Software\\ZooImage"
		#res <- try(tk2reg.setkey(ZIkey), silent = TRUE)
		#assignTemp("ZIkey", ZIkey)
	#}

	## Load the various image resources
	#if (!redef && interactive()) ImgReadPackage("zooimage")

	## Load the menus
	#if (!redef && interactive()) MenuReadPackage("zooimage")

	## Possibly create the ZIguiPackage variable to indicate from where to load
	## other GUI resources
	ZIguiPackage <- getTemp("ZIguiPackage")
	if (!redef || is.null(ZIguiPackage))
		ZIguiPackage <- "zooimage"
	assignTemp("ZIguiPackage", ZIguiPackage)

	## The directory that contains binary executables
	bindir <- system.file("bin", package = "zooimage")

	## Determine where to find ImageJ
	## TODO... currently, it is in a fixed position
	## TODO: no need to ship the exe file, we can just ship a simple
	## bat file with java -jar ij.jar -ijpath=./plugins
	if (interactive()) {
		if (isWin()) {
			ImageJExe <- file.path(bindir, "ImageJ", "ImageJ.exe")
		} else if (isMac()) {
			#ImageJExe <- "/Applications/Fiji/Fiji.app/Contents/MacOS/fiji-macosx"
			ImageJExe <- "open /Applications/Fiji/Fiji.app"
		} else {
			## TODO... Get ImageJ executable
			ImageJExe <- "fiji"
		}
		if (file.exists(ImageJExe)) options(ImageEditor = ImageJExe)
	} else options(ImageEditor = "")

	## Determine where to find XnView
	## TODO... currently, it is in a fixed position
	if (interactive()) {
		if (isWin()) {
			XnViewExe <- file.path(bindir, "XnView", "XnView.exe")
		} else if (isMac()) {
			XnViewExe <- "/Applications/Utilities/XnViewMP.app/Contents/MacOS/xnview"
		} else {
			XnViewExe <- "nautilus --geometry 600x600"
		}
		if (file.exists(XnViewExe)) options(ImageViewer = XnViewExe)
	} else options(ImageViewer = "")
	
	## Determine where to find VueScan
	## TODO... currently, it is in a fixed position
	if (interactive()) {
		if (isWin()) {
			VueScanExe <- file.path(bindir, "VueScan", "VueScan.exe")
		} else if (isMac()) {
			VueScanExe <- "/Applications/VueScan.app/Contents/MacOS/VueScan"
		} else {
			## TODO: other locations for Mac or Linux?!
			VueScanExe <- "vuescan"
		}
		if (file.exists(VueScanExe)) options(VueScan = VueScanExe)
	} else options(VueScan = "")
	
	## Possibly load the ZooImage assistant
	LoadIt <- getOption("ZIAssistant")
	if (is.null(LoadIt) || LoadIt == TRUE) ZIDlg()

	## Set the default template directory
	if (is.null(getOption("ZITemplates")))
		options(ZITemplates = system.file("templates", package = "zooimage"))

	## Switch to the default directory, if defined
	defdir <- getOption("ZI.DefaultDirectory", "")
    if (defdir != "" && file.exists(defdir) && file.info(defdir)$isdir)
        setwd(defdir)
}

## Unloading ZooImage
.onUnload <- function (libpath)
{
	## Eliminate the ZooImage menu entries
	if (.Platform$GUI[1] == "Rgui") {
		try(menuDel("$ConsoleMain/ZooImage"), silent = TRUE)
		try(menuDel("$ConsolePopup/ZooImage"), silent = TRUE)
	}
	closeAssistant()
}

## R version < 2.15.0 does not have paste0 => create it here
if (compareRVersion("2.15.0") < 0) {
	paste0 <- function (..., collapse = NULL)
		paste(..., sep = "", collapse = collapse)
}