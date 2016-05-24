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

## Transforms a file extension to a pattern for ignore.case matching of the  
## extension: extension (with or without the dot at the beginning)
## returns a regular expression pattern that can be used
## to match files with this extension
extensionPattern <- function (extension = "r",
add.dot = !grepl("[.]", extension))
{
	extensionLetters <- substring(extension, 1:nchar(extension),
		1:nchar(extension))
	parts <- ifelse(extensionLetters %in% c(letters, LETTERS), 
		paste("[", extensionLetters, casefold(extensionLetters, upper = TRUE),
		"]", sep = ""), extensionLetters)
	pattern <- paste(parts, collapse = "") 
	if (isTRUE(as.logical(add.dot)))
		pattern <- paste(".", pattern, sep = "")
	pattern <- gsub( "[.]", "[.]", pattern)
	paste(pattern, "$", sep = "")
}

## Checks if the file has the given extension (used at different places...)
hasExtension <- function (file, extension = "r",
pattern = extensionPattern(extension))
	grepl(pattern, file)

## Get the name of a file, without its extension
noExtension <- function (file)
	sub("\\.[^.]+$", "", basename(file))

## List files with given extension
listFilesExt <- function (dir, extension = "r",
pattern = extensionPattern(extension), ... )
{
	if (!checkDirExists(dir)) return(character(0))
	list.files(dir, pattern = pattern , ...)
}

zimList <- function (dir, ...)
	listFilesExt(dir, extension = "zim", ...)

zimDatList <- function (dir, ...)
	listFilesExt(dir, extension = "_dat1.zim", ...)

zipList <- function (dir, ...)
	listFilesExt(dir, extension = "zip", ...)

zidList <- function (dir, ...)
	listFilesExt(dir, extension = "zid", ...)
	
zidbList <- function (dir, ...)
	listFilesExt(dir, extension = "zidb", ...)

jpgList <- function (dir, ...)
	listFilesExt(dir, extension = "jpg", ...)
	
pngList <- function (dir, ...)
	listFilesExt(dir, extension = "png", ...)

## Check if a file exists
checkFileExists <- function (file, extension, message = "file not found: %s",
force.file = FALSE)
{
	## Does this file exists?
	if (!all(file.exists(file))) {
		warning(sprintf(message, file))
		return(FALSE)
	}
	
	## Make sure it is not a directory
	if (force.file && any(file.info(file)$isdir)) {
		warning("one or more files are directories")
		return(FALSE)	
	}
	
	## Check its extension
	if (!missing(extension) && !all(hasExtension(file, extension))) {
		warning(sprintf("one or more files are not '%s' file", extension))
		return(FALSE)
	}
	
	## Everything is fine!
	return(TRUE)
}

## Checks if a directory exists
checkDirExists <- function (dir,
message = 'Path "%s" does not exist or is not a directory')
{
	if (!all(file.exists(dir)) || !all(file.info(dir)$isdir)) {
		warning(sprintf(message, dir))
		return(FALSE)
	}
	
	## Everything is fine...
	return(TRUE)
}

## Check if a directory is empty (used in prepareTrain())
checkEmptyDir <- function (dir, message = 'dir "%s" is not empty')
{	
	## Works only on a single dir (not vectorized code)
	dir <- as.character(dir)[1]
	if (file.exists(dir)) {
		Files <- list.files(dir, all.files = TRUE)
		Files <- Files[!Files %in% c(".", "..")]
		if (length(Files > 0)) {
			warning(sprintf(message, dir))
			return(FALSE)
		} else return(TRUE)
	} else forceDirCreate(dir)	
}

## Force creation of a directory
forceDirCreate <- function (dir)
{	
	## If it exists, make sure it is a directory
	if (file.exists(dir)) {
		if (!file.info(dir)$isdir) {
			warning(sprintf('"%s" is not a directory', dir))
			FALSE
		} else TRUE
	} else if (!dir.create(dir, showWarnings = FALSE)) {
		warning(sprintf('could not create directory "%s"', dir))
		FALSE
	} else TRUE
}

## Checks the first line of a file against some expected content
checkFirstLine <- function (file, expected = c("ZI1", "ZI2", "ZI3"), 
message = 'file "%s" is not a valid ZooImage version <= 3 file')
{
	Line1 <- scan(as.character(file)[1], character(), nmax = 1, quiet = TRUE)
	res <- Line1 %in% expected
	if (!res) warning(sprintf(message, file))
	return(res) 
}
