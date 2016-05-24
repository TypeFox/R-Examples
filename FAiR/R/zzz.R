#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

.onLoad <- function(lib, pkg) {
	library.dynam("FAiR", pkg, lib)
}

.onUnload <- function(libpath) {
	library.dynam.unload("FAiR", libpath)
}

.onAttach <- function( ... ) {
FAiRLib <- dirname(system.file(package = "FAiR"))
version <- packageDescription("FAiR", lib.loc = FAiRLib)$Version
BuildDate <- packageDescription("FAiR", lib.loc = FAiRLib)$Date
packageStartupMessage(paste("##  FAiR Version", version, "Build Date:", BuildDate))
packageStartupMessage("## See http://wiki.r-project.org/rwiki/doku.php?id=packages:cran:fair for more info")
packageStartupMessage("FAiR  Copyright (C) 2008 -- 2012  Benjamin King Goodrich")
packageStartupMessage("This program comes with ABSOLUTELY NO WARRANTY.")
packageStartupMessage("This is free software, and you are welcome to redistribute it")
packageStartupMessage("under certain conditions, namely those specified in the LICENSE file")
packageStartupMessage("in the root directory of the source code.")

if (length(grep("darwin", R.version$platform))) {
  packageStartupMessage("\n\nWARNING: It appears you are using a Mac.\n",
	"FAiR will CRASH under normal usage unless the X server is running.\n",
	"If you are using the R GUI, click the X icon in the top middle ",
	"region of the GUI and then restart the R GUI.\n",
	"If you are not using the R GUI, close R, execute 'open -a X11.app' in ",
	"a shell and then restart R\n",
	"If you do not have an X server installed, it is very difficult to use ",
	"FAiR on a Mac; see the OSX installation disc to install it if necessary.")
}

else if (.Platform$OS.type == "windows") {
	flush.console()
	packageStartupMessage("\n\n It appears you are using Windows.\n",
	    "It is recommended that you disable buffering by pressing Ctrl-W or\n",
	    "by deselecting Misc -> Buffered output in the menu at the top.\n",
	    "Doing so will consistently print the progress of the genetic algorithm.\n")
	flush.console()
}

invisible(NULL)
}

