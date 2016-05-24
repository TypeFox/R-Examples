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

## Zip a .tif image and embed the corresponding .zim file as comment
## This requires the 'zip' program!
zipImg <- function (imagefile, zimfile = NULL, check.zim = TRUE,
replace = FALSE, delete.source = FALSE)
{
	## We need to switch to the image dir for correct path in the zip file
	imagefile <- as.character(imagefile)
	if (length(imagefile) != 1) {
		warning("you must provide exactly one image file name")
		return(invisible(FALSE))
	}
	## Check if imagefile exists
	if (!checkFileExists(imagefile, force.file = TRUE,
		message = "%s doesn't exist, or is a directory!"))
		return(invisible(FALSE))
		
	## Switch directory to the one of the image
	initdir <- setwd(dirname(normalizePath(imagefile)))
	on.exit(setwd(initdir))
	## Simplify image file path, since we are now in the right dir
	imagefile <- basename(imagefile)

	## Is there an associated .zim file?
	if (!length(zimfile)) {
		fraction <- sampleInfo(imagefile, "fraction",
			ext = extensionPattern("tif"))
		zimfile <- paste(fraction, "zim", sep = ".")
	} else {
		zimfile <- as.character(zimfile)
		if (length(zimfile) > 1) {
			warning("you cannot provide more than one ZIM file")
			return(invisible(FALSE))
		}
	}

	### TODO: the zim file can be other parts of it , like Sample+A1.zim,
	###       instead of Sample+A.zim!
	if (!checkFileExists(zimfile, force.file = TRUE,
		message = "%s doesn't exist; cannot process the corresponding image"))
		return(invisible(FALSE))

	## Verify the content of the .zim file (returns -1 in case of error)
	if (isTRUE(as.logical(check.zim)) && zimVerify(zimfile) < 0)
		return(invisible(FALSE))

	## Zip the image in the '_raw' subdir and add the information from the .zim
	## file as comment
	message("Zipping image '", imagefile, "' ...")
	zipfile <- paste(noExtension(imagefile), "zip", sep = ".")
	zipfile <- file.path(".", "_raw", zipfile)
	## Make sure that "_raw" subdir exists
	if (!forceDirCreate("_raw")) return(invisible(FALSE))

	## Copy or move the image to a .zip compressed file
	if (isTRUE(as.logical(replace)) && file.exists(zipfile))
		unlink(zipfile)
	
	## zip() function returns status zero if everything is fine
	if (zip(zipfile, imagefile, flags = "-rq9X") != 0) {
		warning("error while zipping '", basename(imagefile), "'")
		return(invisible(FALSE))
	}

	## Add comment to the zip file
	## Note: the .zim file is never deleted, because it can be used for other
	## purposes!
	## Note2: except for a warning, we don't care about not adding .zim data
	if (!zipNoteAdd(zipfile, zimfile)) {}

	## Do we delete source image? (not much a problem if it fails too)
	if (isTRUE(as.logical(delete.source))) unlink(imagefile)

	## Invisibly indicate success
	invisible(TRUE)
}

## Compress all .tif images in the corresponding directory
## (at least those with an associated .zim file)
zipImgAll <- function (path = ".", images = NULL, check.zim = TRUE,
replace = FALSE, delete.source = FALSE)
{
	## First, switch to that directory
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- "."	# Indicate we are now in the right path

	## Get the list of images to process
	if (!length(images))	# Compute them from path
		images <- dir(path, pattern = extensionPattern("tif")) # All .tif files

	## If there is no images in this dir, exit now
	if (!length(images)) {
		warning("There is no images to process in ", getwd())
		return(invisible(FALSE))	
	}

	## Make sure there is no path associated
	if (!all(images == basename(images))) {
		warning("You cannot provide paths for 'images', just file names")
		return(invisible(FALSE))
	}

	## Look at associated .zim files
	zimfiles <- paste(sampleInfo(images, "fraction",
		ext = extensionPattern("tif") ), ".zim", sep = "")
	keep <- file.exists(zimfiles)
	if (!any(keep)) {
		warning("You must create ZIM files first (ZooImage Metadata)!")
		return(invisible(FALSE))	
	}
	if (!all(keep)) {
    	warning(sum(!keep), " on ", length(keep),
			" images have no ZIM file associated and will not be processed!")
		images <- images[keep]
		zimfiles <- zimfiles[keep]
	}

	## Check the zim files
	ok <- TRUE
	if (isTRUE(as.logical(check.zim))) {
		message("Verification of ZIM files...")
		flush.console()
		zfiles <- unique(zimfiles)
		zimCheck <- function (zim) {
			message("Verifying '", basename(zim), "' ...")
			zimVerify(zim) >= 0
		}
		ok <- batch(zfiles, zimCheck, verbose = FALSE)
	}
	if (ok) {
		message("-- Done! --")
	} else {
		warning("corrupted ZIM file(s) found, compression not started!")
		return(invisible(FALSE))
	}

	## If everything is ok compress these files
	message("Compression of images...")
	flush.console()
	ok <- batch(images, zipImg, check.zim = FALSE,
		replace = replace, delete.source = delete.source, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(images),
			" images were compressed (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Uncompress .tif image and .zim file from a .zip archive file
unzipImg <- function (zipfile, replace = FALSE, delete.source = FALSE)
{
	# Extract .zim file, .tif file or both from a .zip archive
	zipfile <- as.character(zipfile)
	if (length(zipfile) != 1) {
		warning("you must provide one file path in 'zipfile'")
		return(invisible(FALSE))
	}
	
	## Check if zipfile exists
	if (!checkFileExists(zipfile, force.file = TRUE,
		message = "%s doesn't exist, or is a directory!"))
		return(invisible(FALSE))
	
	## Special case: if dir is _raw, then extract into parent dir (..)
	isRawDir <- basename(dirname(normalizePath(zipfile))) == "_raw"
	
	## Switch directory to the one of the zip archive
	initdir <- setwd(dirname(normalizePath(zipfile)))
	on.exit(setwd(initdir))
	## Simplify zip file path, since we are now in the right dir
	zipfile <- basename(zipfile)
	
	## Determine the name of the corresponding .zim file
	fraction <- sampleInfo(zipfile, "fraction",
		ext = extensionPattern("zip"))
	zimfile <- paste(fraction, "zim", sep = ".")
	if (isRawDir) zimfile <- file.path("..", zimfile)
	
	message("Unzipping '", zipfile, "' ...")
	
	## Do we replace existing .zim files?
	replace <- isTRUE(as.logical(replace))
	if (replace || !file.exists(zimfile)) {
		## Extract data from the zimfile
		if (!length(zipNoteGet(zipfile, zimfile)))
			return(invisible(FALSE))
	}
	
	## Unzip the .tif image
	if (isRawDir) exdir <- ".." else exdir <- "."
	if (!length(tryCatch(unzip(zipfile, overwrite = replace, junkpaths = TRUE,
		exdir = exdir), error = function (e) warning(e),
			warning = function (w) return()))) {
		message("    ... not done!")
		return(invisible(FALSE))
	}

	## Do we delete zip archive? (not much a problem if it fails here)
	if (isTRUE(as.logical(delete.source))) unlink(zipfile)

	## Invisibly indicate success
	invisible(TRUE)
}

## Extract all .zim, .tif or both from .zip files
unzipImgAll <- function (path = ".", zipfiles = NULL, replace = FALSE,
delete.source = FALSE)
{
	## First, switch to that directory
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- "."	# Indicate we are now in the right path

	## Get the list of zip archives to process
	if (!length(zipfiles))	# Compute them from path
		zipfiles <- dir(path, pattern = extensionPattern("zip")) # All .zip

	## If there is no .zip files in this dir, exit now
	if (!length(zipfiles)) {
		warning("There is no ZIP archives to process in ", getwd())
		return(invisible(FALSE))	
	}

	## Make sure there is no path associated
	if (!all(zipfiles == basename(zipfiles))) {
		warning("You cannot provide paths for 'zipfiles', just file names")
		return(invisible(FALSE))
	}

	## Uncompress these files
	message("Decompression of .zip archives...")
	flush.console()
	ok <- batch(zipfiles, unzipImg, replace = replace,
		delete.source = delete.source, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(zipfiles),
			" archives were uncompressed (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}
}
