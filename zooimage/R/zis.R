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

## Read data from a .zis file
zisRead <- function (zisfile = "Description.zis", 
expected.sections = c("Description", "Series", "Cruises", "Stations", "Samples"))
{
    if (!checkFileExists(zisfile, extension = "zis", force.file = TRUE))
		return(NULL)
	if (!checkFirstLine(zisfile)) return(NULL)
	rl <- readLines(zisfile,  encoding = "UTF-8")
	if (!length(rl) > 1) {
		warning("the file is empty or corrupted!")
		return(NULL)
	}
	positions <- grep("^[[].*[]]", rl)
	sections <- sub("^[[](.*)[]]", "\\1", rl[positions])
	if (!all(expected.sections %in% sections)) {
		warning("incorrect ZIS file; it does not contain all expected sections")
		return(NULL)
	}
	start <- positions + 1
	end <- c(tail(positions, -1) - 2, length(rl))
	readData <- lapply(1:length(start), function (i) {
		if (sections[i] == "Description") {
			rx <- "^(.*?)=(.*)$"
			txt <- rl[start[i] : end[i]] 
			variables <- sub(rx, "\\1", txt)
			values <- sub(rx, "\\2", txt)
			out <- data.frame(matrix(values, nrow = 1))
			names(out) <- variables
		} else {
			con <- textConnection(rl[start[i] : end[i]])
			on.exit(close(con))
			out <- read.table(con, sep = "\t", header = TRUE, dec = getDec(),
				blank.lines.skip = FALSE)
			names(out)[1] <- sub("^X\\.", "", names(out)[1])
			out <- out[, !grepl("^X\\.[0-9]+", names(out))]
		}
		return(out)
	})
	names(readData) <- sections
	Samples <- readData[["Samples"]]
	Samples$Date <- as.Date(Samples$Date)
	Series <- readData[["Series"]]
	Cruises <- readData[["Cruises"]]
	Cruises$Start <- as.Date(Cruises$Start)
	Cruises$End <- as.Date(Cruises$End)
	Stations <- readData[["Stations"]]
	Stations$Start <- as.Date(Stations$Start)
	Stations$End <- as.Date(Stations$End)
	Description <- readData[["Description"]]
	
	## Combine all this in a data frame + metadata
	structure(Samples, 
		metadata =  list(Desc = Description, Series = Series, Cruises = Cruises,
		Stations = Stations), class = c("ZIDesc", "data.frame"))
}

## Create a .zis file
zisCreate <- function (zisfile, template = NULL,
edit = TRUE, editor = getOption("fileEditor"), wait = FALSE)
{	
	## Use a ui to get the file name
	if (missing(zisfile) || !length(zisfile) || zisfile == "") {
		zisfile <- dlgInput("Give a name for the new ZIS file:",
			title = "ZIS file creation", default = "Description.zis")$res
		if (!length(zisfile)) return(invisible(FALSE))
		if (!hasExtension(zisfile, "zis"))
			zisfile <- paste(zisfile, ".zis", sep = "")
	}
	
    ## If the file already exists, edit current version
	if (file.exists(zisfile))
		if (isTRUE(edit)) {
			return(zisEdit(zisfile, editor = editor, wait = wait))
		} else return(invisible(TRUE))
	
	## Look for the template
	if (is.null(template))
		template <- file.path(getOption("ZITemplates"), "Description.zis")
	if (!checkFileExists(template, "template '%s' not found", extension = "zis"))
		return(invisible(FALSE))
	## Copy the template into the new file
	file.copy(template, zisfile)

	## Possibly edit this new file
	if (isTRUE(edit)) {
		return(zisEdit(zisfile, editor = editor, wait = wait))
	} else return(invisible(TRUE))
}

## Edit a .zis file
zisEdit <- function (zisfile, editor = getOption("fileEditor"), wait = FALSE, ...)
{
    if (missing(zisfile) || !length(zisfile) || zisfile == "") {
		zisfile <- selectFile("Zis")
		if (zisfile == "") return(invisible(FALSE))
	} else if (!checkFileExists(zisfile,
		message = "the file '%s' is not found!", extension = "zis"))
		return(invisible(FALSE))
	fileEdit(zisfile, editor = editor, wait = wait, ...)
}
