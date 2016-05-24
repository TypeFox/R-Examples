################################################################################
# Copyright (C) 2010 by 52 North                                               #
# Initiative for Geospatial Open Source Software GmbH                          #
#                                                                              #
# Contact: Andreas Wytzisk                                                     #
# 52 North Initiative for Geospatial Open Source Software GmbH                 #
# Martin-Luther-King-Weg 24                                                    #
# 48155 Muenster, Germany                                                      #
# info@52north.org                                                             #
#                                                                              #
# This program is free software; you can redistribute and/or modify it under   #
# the terms of the GNU General Public License version 2 as published by the    #
# Free Software Foundation.                                                    #
#                                                                              #
# This program is distributed WITHOUT ANY WARRANTY; even without the implied   #
# WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU #
# General Public License for more details.                                     #
#                                                                              #
# You should have received a copy of the GNU General Public License along with #
# this program (see gpl-2.0.txt). If not, write to the Free Software           #
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or #
# visit the Free Software Foundation web page, http://www.fsf.org.             #
#                                                                              #
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)                          #
# Created: 2010-06-18                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

?package.skeleton

.path = "/home/daniel/Dropbox/2010_SOS4R/workspace/sos4R"
files <- c(
		paste(.path, "R", "Constants.R",  sep = "/"),
		paste(.path, "R", "Class-OWS.R",  sep = "/"),
		paste(.path, "R", "Class-GML.R",  sep = "/"),
		paste(.path, "R", "Class-SWE.R",  sep = "/"),
		paste(.path, "R", "Class-OM.R",  sep = "/"),
		paste(.path, "R", "Class-SA.R",  sep = "/"),
		paste(.path, "R", "Class-OGC.R",  sep = "/"),
		paste(.path, "R", "Class-SOS.R",  sep = "/"),
		paste(.path, "R", "Class-SML.R",  sep = "/"),
		paste(.path, "R", "Generic-methods.R",  sep = "/"),
		paste(.path, "R", "OWS-methods.R",  sep = "/"),
		paste(.path, "R", "OWS-methods-parsing.R",  sep = "/"),
		paste(.path, "R", "SOS-methods-parsing.R",  sep = "/"),
		paste(.path, "R", "OM-methods.R",  sep = "/"),
		paste(.path, "R", "OM-methods-parsing.R",  sep = "/"),
		paste(.path, "R", "SA-methods.R",  sep = "/"),
		paste(.path, "R", "GML-methods.R",  sep = "/"),
		paste(.path, "R", "SWE-methods.R",  sep = "/"),
		paste(.path, "R", "SML-methods.R",  sep = "/"),
		paste(.path, "R", "GML-methods-parsing.R",  sep = "/"),
		paste(.path, "R", "SA-methods-parsing.R",  sep = "/"),
		paste(.path, "R", "SWE-methods-parsing.R",  sep = "/"),
		paste(.path, "R", "OGC-methods.R",  sep = "/"),
		paste(.path, "R", "PrintShowStructureSummmary-methods.R",  sep = "/"),
		paste(.path, "R", "SOS-methods-util.R",  sep = "/"),
		paste(.path, "R", "SOS-methods.R",  sep = "/"),
		paste(.path, "R", "Defaults.R",  sep = "/")
	)

#
# create to temp folder, then copied manually
#
package.skeleton(name = "sos4R", namespace = TRUE, code_files = files,
		path = "/tmp/")


#
# if the package is attached to the session, this method can be used to add only
# certain files and create documentation for them
#
library("SoDA")
?packageAdd
?promptAll

#
# alternatively: create documentation files manually
#
?prompt
?promptMethods
?promptClass
?promptPackage

# for print methods
promptMethods(f = "print", filename = "/home/daniel/Dokumente/2010_SOS4R/workspace/sos4R/man/print-methods.Rd")

# 
promptClass(clName = "OmObservationCollection",
		filename = "/home/daniel/Dokumente/2010_SOS4R/workspace/sos4R/man/OmObservationCollection.Rd")

#
promptMethods(f = "sosFeatureIds", filename = "/tmp/sosFeatureIds-methods.Rd")
promptMethods(f = "sosParse", filename = "/tmp/sosParse-methods.Rd")



#
#
#
library("tools")
?codoc # Check Code/Documentation Consistency

################################################################################
# tools::showNonASCII(readLines('sos4R.Rnw')) 
checkNonASCII <- function(pkgPath) {
	require("tools")
	
	# get all files in the workspace
	p <- pkgPath
	ps <- c("", "demo", "inst", "inst/doc", "man", "R", "sandbox", "tests")
	dirs <- paste(p, ps, sep = "/")
	
	?showNonASCII
	?readLines
	?dir
	
	filenames <- lapply(dirs, dir)
	
	filepaths <- list()
	for (i in seq(along = dirs)) {
		.l <- paste(dirs[[i]], filenames[[i]], sep = "/")
		filepaths <- c(filepaths, .l)
	}
	
	# remove some folders
	filepaths <- filepaths[!grepl(pattern = "//", x = filepaths)]
	filepaths <- filepaths[!grepl("/R$", filepaths)]
	filepaths <- filepaths[!grepl("/inst/doc$", filepaths)]
	filepaths <- filepaths[!grepl(".RData$", filepaths)]
	filepaths <- filepaths[!grepl(".pdf$", filepaths)]
	filepaths
	
	# check characters
	for (i in seq(along = filepaths)) {
		cat(filepaths[[i]], "\n")
		.file <- readLines(filepaths[[i]])
		showNonASCII(.file)
	}
}

checkNonASCII("D:/workspace/sos4R")

################################################################################
# tools::compactPDF
#?tools::compactPDF

# run this before every commit...
result <- tools::compactPDF(paths = "D:/workspace/sos4R/inst/doc")
result
# or even better: run R CMB build with option "--compact-vignettes"
