## Copyright (c) 2006-2012, Ph. Grosjean <phgrosjean@sciviews.org>
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

## Make .zim files and import images, using a .zie import file for specifs
zieMake <- function (path = ".", Filemap = "Import_Table.zie", check = TRUE,
replace = FALSE, move.to.raw = TRUE, zip.images = "[.]tif$")
{
	## Example of use:
	## Import Digicam RAW files (currently, only Canon .CR2 files)
	## and transform them into .pgm file with correct names in _work subdirectory
	## move processed .cr2 files into _raw; create associated .zim files

	## This requires the 'dc_raw' and 'ppmtopgm' programs plus a couple of others!
	## TODO: change this to eliminate external programs dependencies!
	## We need 'identify' and 'convert' from ImageMagick 16 bits!
	## Make sure they are available
	if (isTRUE(as.logical(check))) {
		#checkCapable("identify")
		#checkCapable("convert")
		#checkCapable("dc_raw")
		#checkCapable("ppmtopgm")
		#checkCapable("zip")
	}
	
	## First, switch to the root directory
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- getwd()	# Indicate we are now in the right path
	### TODO If last subdir of path is "_raw", then, work with parent dir
	## and do not move files in _raw subdir
	
	## Read the Filemap
	cat("Reading Filemap...\n")
	if (!checkFileExists(Filemap, extension = "zie", force.file = TRUE))
		return(invisible(FALSE))
	
	## Check first line for ZI1-3
	if (!checkFirstLine(Filemap)) return(invisible(FALSE))
	
	## Read the file and check it is not empty
	## Note: we don't use comment.char = '#' because we want to read and rewrite
	## those comments!
	Lines <- scan(Filemap, character(), sep = "\t", skip = 1,
		blank.lines.skip = FALSE, flush = TRUE, quiet = TRUE, comment.char = "") 
	if (!length(Lines)) {
		warning('filemap empty or corrupted!')
		return(invisible(FALSE))
	}
	
	## Get the position of a section
	getSectionPos <- function (section = "Map",
		message = "section '[%s]' found") {
		rx <- sprintf("[[]%s[]]", section)
		out <- grep(rx, Lines)
		if (length(out) != 1) {
			warning(sprintf(message, section))
			NULL
		} else out
	}
	
	getSection <- function (section = "Map", to = c("next","end"),
		message = "The [Map] section is empty!") {
		to <- match.arg(to)
		start <- getSectionPos(section)[1]
		if (!length(start)) return(NULL)
		end <- switch(to, 
			"next" = {
				ends <- getSectionPos(".*")
				if (!length(ends)) return(NULL)
				ends[ends > start][1] - 1
			}, 
			"end" = length(Lines)
		)
		out <- Lines[seq.int(from = start + 1, to = end)]
		if (!length(out)) {
			warning(message)
			NULL
		} else out
	}
	
	## Get everything before '[Map]' as template data for the .zim file
	posMap <- getSectionPos("Map",
		"The file is corrupted: no or duplicated [Map] section found!")
	if (!length(posMap)) return(invisible(FALSE))
		
	## Setup the zim data
	zimData <- Lines[1:(posMap - 1)]
	attr(zimData, "Sample") <- NULL	# Currently, there is no sample!
	attr(zimData, "MakeZim") <- FALSE
	
	## Extract various properties 
	
	## Property extractor
	property <- function (property = "FilenamePattern", default = "") {
		rx <- sprintf("^%s[[:space:]]*=[[:space:]]*(.*)", property)
		if (any(gl <- grepl(rx, Lines))) {
			sub(rx, "\\1", Lines[gl][1])
		} else default
	}
	
	FilePat <- property("FilenamePattern")
	FracPat <- property("FractionPattern")
	SubPat <- property("SubsamplePattern")
	Convert <- property("Convert")
	Return <- property("Return")
	FileExt <- property("FileExt")
	FileC <- property("FileC")
	FileExt2 <- property("FileExt2", FileExt)
	MoveToWork <- tolower(property("MoveToWork")) %in% c("true", "yes", "1")
	Exif <- property("[<]exif[>]") != ""
	attr(zimData, "Exif") <- "" # Nothing yet here
	
	## Get the [Map] section
	Lines <- getSection("Map", to = "end", "The [Map] section is empty!")
	if (!length(Lines)) return(invisible(FALSE))
	
	message("Reading Filemap... OK!")
	
	## Make sure _raw, and _work subdirectories exists and have write access
	if (!forceDirCreate("_raw")) return(invisible(FALSE))
	if (Convert != "" || MoveToWork)
		if (!forceDirCreate("_work")) return(invisible(FALSE))
	
	## This function constructs image filename using possibly a FilenamePattern
	MakeImageName <- function(x, pattern = FilePat) {
		if (pattern == "") return(x)
		
		## Do we need to format a number?
		Format <- sub("^.*[<]([1-9]?)[>].*$", "\\1", pattern)
		if (Format != "")
			x <- formatC(as.integer(x), width = as.integer(Format), flag = "0")
		
		## Make the replacement according to the pattern
		File <- gsub("[<][1-9]?[>]", x, pattern) # Do we have to use FilePattern?
		return(File)
	}
	
	## Make sure that all image files are there, and there is no duplicated use
	## of the same image
	### TODO: indicate progression with exact line number in the zie file!
	### TODO: allow restarting from a given point!
	message("Checking all lines in the ZIE file for raw images...")
	allImages <- character(0)
	nLines <- length(Lines)
	for (i in 1:nLines) {
		### TODO: allow restarting from a given point and eliminate previous 
		###       lines for which there are no images (considered as already
		###       processed!)
		progress(i, nLines)
		if (!grepl("^[-][>]", Lines[i])) {	# This is not a state change command
			File <- MakeImageName(trimString(sub("[=].*$", "", Lines[i])))
			checkFileExists(File)
			if (File %in% allImages) {
				warning(sprintf("Duplicated use of the same file : '%s' !",
					File))
				return(invisible(FALSE))
			}
			allImages <- c(allImages, File)		
		}
	}
	progress(101) # Clear progression indicator
	message("...OK!")
	
	## Now that we know all image files are there, process the [Map] section
	## line-by-line	
	message("Processing all lines in the ZIE file (import images and make ZIM files)...")
	ok <- TRUE
	
	## BuildZim : This function builds the zim file and check it
	BuildZim <- function (zimData, FracPat, SubPat) {	
		## Calculate the name of the zim file
		zimFileName <- paste(Smp, "zim", sep = ".")
		zimFile <- file.path(getwd(), zimFileName)
		
		## If the zim file already exists, skip this
		if (!replace && file.exists(zimFile)) {
			warning(".zim file already exists for '", Smp, "'")
			return(TRUE)
		}
		
		## Make necessary replacement in Fraction and Subsample
		Smp <- attr(zimData, "Sample")
		if (is.null(Smp) || Smp == "") return(FALSE)
		zimD <- zimData
		
		## Clear a whole section, starting from its header to the next header
		ClearSection <- function (Data, fromLine) {
			n <- length(Data)
			if (fromLine > n) return(Data)
			## Locate the next header (line starting with "[")
			NextHeader <- grep("^[[]", Data[(fromLine + 1):n])
			if (length(NextHeader) == 0) {
				toLine <- n
			} else {
				toLine <- NextHeader[1] + fromLine - 1
			}
			## Strip out this section
			return(Data[-(fromLine:toLine)])
		}
		
		## Deal with FracPat
		if (FracPat != "") {
			## This is the header to consider
			if (length(grep(FracPat, Smp)) == 0) {
				stop( paste("Sample '", Smp,
					"' is incompatible\nwith FractionPattern '", FracPat, "'",
					sep = ""))
			}				
			Frac <- paste("[[]Fraction_", sub(FracPat, "\\1", Smp), "\\]",
				sep = "")
			posFrac <- grep(Frac, zimD)
			if (length(posFrac) < 1) {
				warning("[Fraction] section not found (", Frac, ")!")
				return(FALSE) 
			}			
			if (length(posFrac) > 1) {
				warning("multiple", Frac, "sections for sample '", Smp, "'")
				return(FALSE)
			}
			zimD[posFrac] <- "[Fraction]"
			## Strip out all other [Fraction_XXX] sections
			otherFrac <- grep("[[]Fraction_", zimD)
			if (length(otherFrac) > 0) 
				for (i in 1:length(otherFrac))
					zimD <- ClearSection(zimD, otherFrac[i])
		}
		
		if (SubPat != "") {
			## This is the header to consider
			if (!length(grep(SubPat, Smp))) {
				warning("Sample '", Smp,
					"' is incompatible\nwith SubsamplePattern '", SubPat, "'")
				return(FALSE)
			}
			Sub <- paste("[[]Subsample_", sub(SubPat, "\\1", Smp), "\\]",
				sep = "")
			posSub <- grep(Sub, zimD)
			if (!length(posSub)) {
				warning("[Subsample] section not found (", Sub, ")!")
				return(FALSE)
			}
			if (length(posSub) > 1) {
				warning("multiple", Sub, "sections found for this sample!")
				return(FALSE)
			}
			zimD[posSub] <- "[Subsample]"
			## Strip out all other [Subsample_XXX] sections
			otherSub <- grep("[[]Subsample_", zimD)
			if (length(otherSub) > 0) 
				for (i in 1:length(otherSub))
					zimD <- ClearSection(zimD, otherSub[i])
		}
		## Possibly insert Exif data
		if (Exif && !is.null(attr(zimData, "Exif"))) {
		    pos <- grep("^[<]exif[>]", zimD)
			# pos is recalculated here, because it may have changed!
		    if (length(pos) > 0)
				zimD <- c(zimD[1:(pos - 1)], attr(zimData, "Exif"),
					zimD[(pos+1):length(zimD)])
		}
				
		## Write the zim file
		message("Writing .zim file for sample '", Smp, "'")
		cat(paste(c("ZI1", zimD), collapse = "\n"), file = zimFile)
		return(TRUE)
	}
	
	## UpdateZim ; This function looks if the line asks for updating zimData and
	## does it (returns TRUE), or it returns FALSE
	UpdateZim <- function (dat, zimData) {
		### TODO: Strip out comments (not done here, because we want to process
		### strings with '#' correctly!
		if (length(grep("^[-][>]", dat)) == 0) return(NULL)
		## This line starts with "->" => we update zimData
		Key <- sub("^[-][>]([^ =]+).*$", "\\1", dat)
		## Special treatment if Key == "Sample"
		if (Key == "Sample") {
			attr(zimData, "Sample") <- trimString(sub("^[^=]+=", "", dat))
			## Indicate that we process another sample
			attr(zimData, "MakeZim") <- TRUE # Tell to make the zim file
			attr(zimData, "Exif") <- ""
		} else { # This is an usual key
			## Replace every line corresponding to this key in zimData
			MatchLines <- grep(paste("^", Key, sep = ""), zimData)
			if (length(MatchLines > 0))
				zimData[MatchLines] <- sub("^[-][>]", "", dat)
		}
		return(zimData)	
	}
	
	## SetCalib : Add or change an entry in [Calibration] section
	SetCalib <- function (Data, Key, Entry) {
		Line <- paste(Key, Entry, sep = "=")
		## Is this key already defined?
		posKey <- grep(paste("^\\s*", Key, "\\s*=", sep = ""), Data)
		## If defined => change it now
		if (length(posKey) > 0) {
			Data[posKey] <- Line
			return(Data)
		}
		## Is the [Calibration] section already defined?
		posCalib <- grep("[[]Calibration\\]", Data)
		if (length(posCalib) > 0) {
			## Add this new entry in the [Calibration] section
			Data <- c(Data[1:posCalib[1]], Line,
				Data[(posCalib[1] + 1):length(Data)])
		} else {
			## Create the [Calibration] section at the end and add this entry
			## inside it
			if (Data[length(Data)] != "")
				Data <- c(Data, "")
			## Make sure that the section is separated with a blank
			Data <- c(Data, "[Calibration]", Line)
		}
		return(Data)
	}
	   
	## Main Loop
	BlankField <- NULL  # The name of the blank-field image to use
	for (i in 1:nLines) {
		progress(i, nLines)
		res <- UpdateZim(Lines[i], zimData)
		if (!length(res)) {
			warning("problem while updating .zim files")
			return(invisible(FALSE))
		}
		
		## This is not a state change command
		if (length(res) == 1 && res == FALSE) {	
			File <- MakeImageName(trimString(sub("[=].*$", "", Lines[i])))
			
			## Determine the name of the converted file
			if (Convert != "") {
				if (FileC == "") { # Construct the name of the converted file
					FileConv <- paste(noExtension(File), FileExt, sep = ".")
				} else {
					## Make sure that previous file is deleted
					unlink(FileC)
					FileConv <- FileC
				}
			} else {
				FileConv <- File
			}
			
			## Determine the final name to give to this converted file,
			## and check if it is a calibration file
			FileConvExt <- tolower(sub("^.*[.]", "", FileConv))
			## Calculate the final name we want for the converted file
			NewFile <- trimString(sub("^.*[=]", "", Lines[i]))
			## 1) If this is 'key' or 'key=' (NeWFile == ""), then,
			##    the file is not renamed!
			if (NewFile == "") {
				FileConvName <- paste(noExtension(File), FileExt2, sep = ".")
				## 2) If the new name starts with "_Calib", then, never use the
				##    Sample part and add a CalibXX entry in .zim file
			} else if (length(grep("^_Calib", NewFile)) > 0) {
                ## If this is a blank-field image, use it for further process
				if (length(grep("^_CalibBF", NewFile)) > 0) {
					FileConvName <- paste(NewFile, FileExt, sep = ".")
					## Remove previous blank-field from root directory
					## (not needed any more!)
					if (!is.null(BlankField)) {
						## Delete blank-field images (.pgm and .img)
						## in the root directory
						unlink(BlankField)
						unlink(paste(noExtension(BlankField), "img", sep = "."))
					}
					BlankField <- FileConvName
				} else {
				    FileConvName <- paste(NewFile, FileExt2, sep = ".")
				}
				## Add or change the calibration information
				## (_CalibSP01 => CalibSP=_CalibSP01.ext)
				Key <- sub("^_(Calib[A-Z]+).*$", "\\1", NewFile)
				zimData <- SetCalib(zimData, Key, FileConvName)
				## 3) Name is Sample + name + ext
			} else {
				Smp <- attr(zimData, "Sample")
				if (is.null(Smp)) Smp <- ""
				FileConvName <- paste(Smp, NewFile, ".", FileExt2, sep = "")
			} 

			## Possibly read Exif data and place it in the zim file
			## (or check correspondance)
			if (Exif) {
				ExifData <- attr(zimData, "Exif")
				ExifData2 <- .readExifRaw(File, check = FALSE)
				if (!is.null(ExifData) && length(ExifData) > 0 &&
					ExifData != "") { # Do a comparison of Exif data
				    compa <- .compareExif(ExifData, ExifData2)
				    if (length(compa) > 0)
						warning("Exif seems to be different from the rest in '",
							File, "'")
				} else { # Just set Exif data
				    attr(zimData, "Exif") <- ExifData2
				}
			}
			
			## Possibly write a zim file?
			MakeZim <- attr(zimData, "MakeZim")
			if (!is.null(MakeZim) && MakeZim) {
				if (BuildZim(zimData, FracPat, SubPat)) {
					attr(zimData, "MakeZim") <- FALSE
				} else {
					return(invisible(FALSE))		
				}
			}
			
			## Possibly convert this file
			if (Convert != "") {
				if (zip.images != "" &&
					length(grep(zip.images, FileConvName)) != 0 &&
					length(grep("^_Calib", FileConvName)) == 0) {
					finalname <- paste(noExtension(FileConvName), "zip",
						sep = ".")
				} else finalname <- FileConvName
				message("Converting image '", File, "' into '", finalname, "'")
				if (replace || !file.exists(FileExt)) { 
					## Create variables Rawbase and Rawnoext
					Rawbase <- File
					Rawnoext <- noExtension(File)
					## Run the command
					res <- eval(parse(text = Convert))
					if (Return != "" && length(grep(Return, res)) == 0) {
						## Check that result matches
						ok <- FALSE
						warning("result after conversion does not match for '",
							File, "'")
					}
					## Look if the converted file is created
					if (!file.exists(FileConv)) {
						ok <- FALSE
						warning("problem: converted file not found '", File, "'")
					}
				}
			} else {
				if (Return != "")
					message("Processing image '", File, "'")
			}

			## If this is a blank-field, then test it
            if (length(grep("^_CalibBF", NewFile)) > 0) {
				msg <- .checkBF(FileConv)
				if (!is.null(msg) && length(msg) > 0 && msg != "") {
					warning(paste(c(
						"Warning! Problem(s) detected with blank-field image:",
						msg), collapse = "\n\t"))  # Report the problem
				}
				## Eliminate dusts and smooth the image with a median filter
				## on a 10 times reduced version of the image
				## We need identify and convert form ImageMagick...16
####				Size <- imagemagick_identify(FileConv)
				Size <- 0  ### TODO: calculate this differently!
				Size2 <- round(Size / 10) # size of the resized image
####				imagemagick_convert(FileConv, Size, Size2)
				
			} else { # make blank-field correction
			    if (!is.null(BlankField)) {
					tryCatch({
						.BFcorrection(FileConv, BlankField, deleteBF = FALSE)
						}, error = function (e) {
							warning(as.character(e))
						})
					
					## Delete the uncorrected file
					unlink(FileConv)
					
					## Now, FileConv is the same file, but with a .tif extension
					FileConv <- paste(noExtension(FileConv), "tif", sep = ".")
					if (!file.exists(FileConv)) {
						ok <- FALSE
						warning("problem: blank-field corrected file not found: '",
							File, "'")
					}
			    }
			}
			
			## If this is an optical density calibration, proceed with it
			if (length(grep("^_CalibOD", NewFile)) > 0) {
				Cal <- calibrate(FileConv)
				Msg <- attr(Cal, "msg")
				## Report the problem
				if (!is.null(Msg) && length(Msg) > 0 && Msg != "") {
					warning(paste(c("problem(s) detected with O.D. calibration image:",
						attr(Cal, "msg")), collapse = "\n\t"))  
				}
				## Put calibration data in the .zim file
				zimData <- SetCalib(zimData, "WhitePoint", round(Cal[1]))
                zimData <- SetCalib(zimData, "BlackPoint", round(Cal[2]))
			}
			
			### TODO: do the same for the spatial calibration image...
			if (Convert == "") {
				## If a second extention is provided, we need to rename and
				## place the original into _raw subdir
				if (FileExt2 != "" && FileExt2 != FileExt) {
					## Copy the original image (indeed, same image, but with
					## original name) into _raw
                	RawFile <- file.path(getwd(), "_raw", File)
                	file.copy(File, RawFile)
                	## And rename the original copy
                	## Possibly move it to _work subdirectory
                	if (MoveToWork)
						FileConvName <- file.path(dirname(FileConvName),
							"_work", basename(FileConvName))
                	file.rename(File, FileConvName)
                }
            } else { # This image was converted
				## Save the original file in _raw subdir
				RawFile <- file.path(getwd(), "_raw", File)
				file.rename(File, RawFile)

				## Rename the converted file and place it in _work
				WorkFileConv <- file.path(getwd(), "_work", FileConvName)
				## Move it, except if it is a blank-field file, then, copy it!
				if (length(grep("^_CalibBF", FileConvName)) > 0) {
					file.copy(FileConv, WorkFileConv)
                	file.rename(FileConv, FileConvName)
				} else {
					file.copy(FileConv, WorkFileConv)
				}
				if (!file.exists(WorkFileConv)) {
					warning("problem moving the converted file into '_work' subdirectory for '",
						File,"'")
					return(invisible(FALSE))
				} else {
					## Do we zip the resulting images, using the zim file
					## as zip comment?
					if (length(grep("^_Calib", FileConvName)) == 0) {
						## Only images, not calib files!
						if (zip.images != "" &&
							length(grep(zip.images, FileConvName)) != 0) {
							curdir <- getwd()
							setwd(file.path(curdir, "_work"))
							zimfile <- paste(attr(zimData, "Sample"), "zim",
								sep = ".")
							# file.copy(file.path(curdir, zimfile), zimfile)
							zipfile <- paste(noExtension(FileConvName), "zip",
								sep = ".")
							zip(zipfile, FileConvName, flags = "-rq9X")
							unlink(FileConvName, recursive = TRUE)
							## Add .zim data as comment to the .zip file
							## Note: except for a warning,
							## we don't care about not adding .zim data
							if (!zipNoteAdd(zipfile,
								file.path(curdir, zimfile))) {}
	
							setwd(curdir)
							## Verify that the .zip file is created
							if (!file.exists(zipfile)) {
								warning(sprintf(
									"problem creating the file : '%s' !",
									zipfile))
								return(invisible(FALSE)) 
							}
						} else {
							### TODO: what do we have to do here???? 
						}
		    		}
				}
			}
		} else zimData <- res
		## Update zimData with value returned by UpdateZim()
	} 
	progress(101) # Clear progression indicator
	
	## Possibly remove latest blank-field from root directory (not needed any more!)
	if (!is.null(BlankField)) {
		## Delete blank-field images (.pgm and .img) in the root directory
		unlink(BlankField)
		unlink(paste(noExtension(BlankField), "img", sep = "."))
	}
	
	if (ok) {
		if (move.to.raw)
			file.rename(Filemap, file.path(getwd(), "_raw", Filemap))
		## There is a bug: a 'fileconv.tif' file is created,
		## delete it for the moment
		unlink("fileconv.tif")
	}
	invisible(TRUE)
}
## example:
## setwd("g:/zooplankton/Madagascar2Macro")	# My example directory
## zieMake(path = ".", Filemap = "Import_Madagascar2Macro.zie")

zieCompile <- function (path = ".", Tablefile = "Table.txt",
Template = "ImportTemplate.zie", Filemap = paste("Import_", noExtension(Tablefile),
".zie", sep = ""), Nrange = c(1, 1000), replace = TRUE, make.it = FALSE)
{		
	message("Creating .zie file...")
	
	## Full path for Filemap
	FilemapPath <- file.path(path, Filemap)
	
	## First, switch to the root directory
	if (!checkDirExists(path)) return(NULL)
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- getwd() # Indicate we are now in the right path
	
    ## Check if needed files exist
    if (!checkFileExists(Tablefile)) return(NULL)
	if (!checkFileExists(Template)) return(NULL)
	
	## Check if the zie file already exists
    if (!isTRUE(as.logical(replace)) && file.exists(Filemap)) {
		warning("'", Filemap,
			"' already exists and is not replaced (replace = FALSE)!")
		return(FilemapPath)
	}
	
	## Read the data from the table
    Data <- read.table(Tablefile, header = TRUE, sep = "\t", dec = getDec(),
		as.is = TRUE)
    
	## Possibly get Nmin and Nmax from the template file
	Nmin <- Nrange[1] # Min number of images for each sample
	Nmax <- Nrange[2] # Max number of images for each sample
    
	## We start from the template
    file.copy(Template, Filemap, overwrite = TRUE)
    Cat <- function (...) cat(..., sep = "", file = Filemap, append = TRUE)
    
	Cat("\n")
    Cat("[Map]\n")
	
	CBF <- -1; CBFNum <- 0
	COD <- -1; CODNum <- 0
	CSp <- -1; CSpNum <- 0
    
	for (i in 1:nrow(Data)) {
        ## Get calibration data
        CalibBF <- Data$CalibBF[i]
        if (!is.na(CalibBF) && !is.null(CalibBF) && CalibBF != "" &&
			CalibBF != CBF) {
			CBFNum <- CBFNum + 1
			text <- paste(CalibBF, sprintf("_CalibBF%3.3d", CBFNum), sep = "=")
			Cat(text, "\n")
			CBF <- CalibBF
        }
        CalibOD <- Data$CalibOD[i]
        if (!is.na(CalibOD) && !is.null(CalibOD) && CalibOD != "" &&
			CalibOD != COD) {
			CODNum <- CODNum + 1
			text <- paste(CalibOD, sprintf("_CalibOD%3.3d", CODNum), sep = "=")
			Cat(text, "\n")
			COD <- CalibOD
        }
		
        CalibSp <- Data$CalibSP[i]
        if (!is.na(CalibSp) && !is.null(CalibSp) && CalibSp != "" &&
			CalibSp != CSp) {
			CSpNum <- CSpNum + 1
			text <- paste(CalibSp, sprintf("_CalibSP%3.3d", CSpNum), sep = "=")
			Cat(text, "\n")
			CSp <- CalibSp
        }
		
        ## Calculate list of all images
        num <- Data$Image[i]
        num <- gsub(";", ",", num, fixed = TRUE)
        num <- gsub("-", ":", num, fixed = TRUE)
        num <- paste("c(", num, ")", sep = "")
        num <- eval(parse(text = num))
        ## Check if the number is correct
		### TODO: add this in the template file!
		if (length(num) < Nmin || length(num) > Nmax) {
			warning("Wrong number of images in 'Image' field for ",
				Data$Sample[i], "!")
			return(NULL)
		}
		
		## Update several fields according to definitions in the samples table
		###TODO: add the other fields + define this option
		Fields <- c("Sample", "SubPart", "PixelSize", "VolIni", "VolPrec")
        Cols <- names(Data)
        for (j in 1:length(Fields)) {
            if (Fields[j] %in% Cols) { 
                value <- Data[i, Fields[j]]
                if (!is.null(value) && !is.na(value) && value != "") {
                    text <- paste("->", Fields[j], "=", value, sep = "")
                    Cat(text, "\n")
                }
            }
        }
		
        ## Insert corresponding images
        for (j in 1:length(num)) {
            text <- paste(num[j], "=.", j, sep = "")
            Cat(text, "\n")
        }
    }
	
	## Do we make it also?
	if (isTRUE(make.it)) {
		res <- zieMake(path = path, Filemap = Filemap, check = TRUE)
		if (res) { # Everything is fine...
			## Move the table and copy the template to the '_raw' subdir too
			file.rename(Tablefile, file.path(path, "_raw", basename(Tablefile)))
			## Move also possibly the .xls equivalent
			Tablexls <- sub("\\.[tT][xX][tT]$", ".xls", Tablefile)
			if (Tablexls != Tablefile && file.exists(Tablexls))
			    file.rename(Tablexls, file.path(path, "_raw",
					basename(Tablexls)))
			file.rename(Template, file.path(path, "_raw", basename(Template)))
		}
	}
	
	## Everything is fine, return the path of the vcreated filemap file
	FilemapPath
}
## example:
## setwd("g:/zooplankton/Madagascar2Macro") # Directory with the example dataset
## zieCompile(Tablefile = "Madagascar2Macro.txt", Nrange = c(2,2))

## Create .zim files and the FitVisParameters.csv file for FlowCAM images
zieCompileFlowCAM <- function (path = ".", Tablefile,
Template = "ImportTemplate.zie", check.names = TRUE)
{
	## Import data from the FlowCAM
	if (!is.character(path) || !file.exists(path) || !file.info(path)$isdir) {
		warning("You must select a path containing text file for FlowCAM images")
		return(invisible(FALSE))
	}
	
	Tablefile <- file.path(path, basename(Tablefile))
	if (!checkFileExists(Tablefile, "txt", force.file = TRUE)) {
		warning("Tablefile not found: '", basename(Tablefile), "'")
		return(invisible(FALSE))
	}
	
	## Read this file
	ImportFile <- read.table(Tablefile, header = TRUE, sep = "\t", dec = ".")

	## Check colnames
	if (isTRUE(as.logical(check.names))) {
		ColNames <- c("Station", "Date", "FlowCell", "Mode", "Magnification",
			"Exp_Name", "Sample", "Dilution", "Sieve", "Volume", "Pump_Speed",
			"Duration", "Temperature", "Salinity", "Gain_Fluo_Ch1",
			"Threshold_Fluo_Ch1", "Gain_Fluo_Ch2", "Threshold_Fluo_Ch2",
			"Threshold_Scatter", "Min", "Max", "Size")
		if (!all(ColNames %in% colnames(ImportFile))) {
			warning("Your import file contains missing columns among Station,",
				" Date, FlowCell, Mode, Magnification, Exp_Name, Sample,",
				" Dilution, Sieve, Volume, Pump_Speed, Duration, Temperature,",
				" Salinity, Gain_Fluo_Ch1, Threshold_Fluo_Ch1, Gain_Fluo_Ch2,",
				" Threshold_Fluo_Ch2, Threshold_Scatter, Min, Max, or Size")
			return(invisible(FALSE))
		}
	}
  
	## Check if the ImportTemplate.zie is present in the directory
	Zie <- file.path(dirname(path), basename(Template))
	if (!file.exists(Zie)) {
		warning("Your directory must contain an 'ImportTemplate.zie' file")
		return(invisible(FALSE))	
	}
  
	## Check if all samples are in the directory and export missing files
	notThere <- character(0)
	for (i in 1:length(ImportFile$Exp_Name))
		if (!file.exists(file.path(dirname(path), ImportFile$Exp_Name[i]))) {
			notThere <- c(notThere, as.character(ImportFile$Exp_Name[i]))
			warning(ImportFile$Exp_Name[i], " is not in the process directory")
		}
  
	## Select only samples present in the process directory
	if (length(notThere)) {
		ImportFile <- ImportFile[!ImportFile$Exp_Name %in% notThere, ]
		warning("import only samples in the process dir or import text file")
	}

	## Read ctx files of the samples from
	Ctx <- file.path(dirname(path), ImportFile$Exp_Name,
		paste(ImportFile$Exp_Name, "ctx", sep = "."))
	lCtx <- length(Ctx)
	if (!lCtx) CtxFile <- NULL else .ctxReadAll(ctxfiles = Ctx)

	## Create fields to generate a table as txt format for the importation
	Experiment <- ImportFile$Exp_Name
	Sample <- ImportFile$Sample
	Image <- CtxFile$Sample_Name
	PixelSize <- CtxFile$pixelsize
	minsize <- CtxFile$minsize
	maxsize <- CtxFile$maxsize
	VolumeDigitized <- CtxFile$VolumeDigitized
	Dilution_VIS <- CtxFile$Dilution_VIS
	SubPart <- ImportFile$Dilution / 100
	VolIni <- ImportFile$Volume
	CellPart <- VolumeDigitized / VolIni
	## Table with value to change
	ImportTxt <- data.frame(Experiment, Sample, Image, PixelSize, SubPart,
		minsize, maxsize, VolumeDigitized, Dilution_VIS, VolIni, CellPart)

	## Read the "ImportTemplate.zie" file
	ZieFileOrig <- scan(Zie, character(), sep = "\t", skip = 0,
		blank.lines.skip = FALSE, flush = TRUE, quiet = TRUE, comment.char = "")
	
	## Loop to create a .zim file
	for (i in 1:nrow(CtxFile)) {
		ZieFile <- ZieFileOrig
		
		## Complete fields using ImportTxt
		ImageLine <- grep("^Sample", ZieFile)
		Sample <- as.numeric(sub("[ ]*$", "", sub("^Sample[ ]*[=][ ]*", "",
			ZieFile[ImageLine[1]])))
		if (is.na(Sample)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$Sample[i], sep = "")
		ImageLine <- grep("^Experiment", ZieFile)
		Experiment <- as.numeric(sub("[ ]*$", "",
			sub("^Experiment[ ]*[=][ ]*", "", ZieFile[ImageLine[1]])))
		if (is.na(Experiment)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$Experiment[i], sep = "")
		ImageLine <- grep("^SubPart", ZieFile)
		SubPart <- as.numeric(sub("[ ]*$", "", sub("^SubPart[ ]*[=][ ]*", "",
			ZieFile[ImageLine[1]])))
		if (is.na(SubPart)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$SubPart[i], sep = "")
		ImageLine <- grep("^CellPart", ZieFile)
		CellPart <- as.numeric(sub("[ ]*$", "", sub("^CellPart[ ]*[=][ ]*", "",
			ZieFile[ImageLine[1]])))
		if (is.na(CellPart)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$CellPart[i], sep = "")
		ImageLine <- grep("^Dilution_VIS", ZieFile)
		Dilution_VIS <- as.numeric(sub("[ ]*$", "",
			sub("^Dilution_VIS[ ]*[=][ ]*", "", ZieFile[ImageLine[1]])))
		if (is.na(Dilution_VIS)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$Dilution_VIS[i], sep = "")
		ImageLine <- grep("^VolIni", ZieFile)
		VolIni <- as.numeric(sub("[ ]*$", "", sub("^VolIni[ ]*[=][ ]*", "",
			ZieFile[ImageLine[1]])))
		if (is.na(VolIni)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$VolIni[i], sep = "")
		ImageLine <- grep("^VolumeDigitized", ZieFile)
		VolumeDigitized <- as.numeric(sub("[ ]*$", "",
			sub("^VolumeDigitized[ ]*[=][ ]*", "", ZieFile[ImageLine[1]])))
		if (is.na(VolumeDigitized)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$VolumeDigitized[i], sep = "")
		ImageLine <- grep("^minsize", ZieFile)
		minsize <- as.numeric(sub("[ ]*$", "", sub("^minsize[ ]*[=][ ]*", "",
			ZieFile[ImageLine[1]])))
		if (is.na(minsize)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$minsize[i], sep = "")
		ImageLine <- grep("^maxsize", ZieFile)
		maxsize <- as.numeric(sub("[ ]*$", "", sub("^maxsize[ ]*[=][ ]*", "",
			ZieFile[ImageLine[1]])))
		if (is.na(maxsize)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$maxsize[i], sep = "")
		ImageLine <- grep("^PixSize", ZieFile)
		PixSize <- as.numeric(sub("[ ]*$", "", sub("^PixSize[ ]*[=][ ]*", "",
			ZieFile[ImageLine[1]])))
		if (is.na(PixSize)) ZieFile[ImageLine[1]] <-
			paste(ZieFile[ImageLine[1]], ImportTxt$PixelSize[i], sep = "")
		
		## Read all context file
		ContextFile <- scan(Ctx[i], character(), sep = "\t", skip = 0,
			blank.lines.skip = FALSE, flush = TRUE, quiet = TRUE,
			comment.char = "")

		## Read note
		Note <- file.path(dirname(path), ImportFile$Exp_Name[i],
			paste(ImportFile$Exp_Name[i], "_notes.txt", sep = ""))
		NoteFile <- scan(Note, character(), sep = "\t", skip = 0,
			blank.lines.skip = FALSE, flush = TRUE, quiet = TRUE,
			comment.char = "")
		
		## Write the resulting table in the sample directory
		Tab <- c(ZieFile, "", ContextFile, "", "[Notes]", NoteFile)
		Export <- file.path(dirname(path), CtxFile$Sample_Name[i],
			paste(CtxFile$Sample_Name[i], "zim", sep = "."))
		write(Tab, file = Export)
	}
	
	## Create a batch file for FlowCAM image analysis using FitVis

	## Select a directory containing a series of FlowCAM runs
	ContextList <- .ctxReadAll(path = path, fill = FALSE, largest = FALSE,
		vignettes = TRUE, scalebar = TRUE, enhance = FALSE, outline = FALSE,
		masks = FALSE, verbose = TRUE)
	
	## Write the table of importation in that directory
	write.table(ContextList, sep = ",", dec = ".", row.names = FALSE, 
		file = file.path(path, "FitVisParameters.csv"), quote = TRUE,
		col.names = TRUE)
	
	message("Import data table has been created in FitVisParameters.csv")

	invisible(TRUE)
}

## The function that eases creation of a ZIE object
### TODO: add a 'message' entry = message to display after of the importation
ZIE <- function (title, filter, description, pattern, command, author,
version, date, license, url, depends = "R (>= 2.4.0), zooimage (>= 1.0-0)",
type = c("import", "export"))
{	
	if (!is.character(title) || !is.character(filter) ||
		!is.character(description) || !is.character(pattern) ||
		!is.character(command) || !is.character(author) ||
		!is.character(version) || !is.character(date) ||
		!is.character(license) || !is.character(url) ||
		!is.character(depends))
		stop("All arguments must be character strings!")
	obj <- list(title = title[1], filter = filter[1], 
		description = paste(description, collapse = "\n"), pattern = pattern[1],
		command = paste(command, collapse = "\n"), author = author[1],
		version = version[1], license = license[1], depends = depends[1])
	type <- match.arg(type, several.ok = FALSE)
	class(obj) <- switch(type,
		import = c("ZIEimport", "ZIE"),
		export = c("ZIEexport", "ZIE"))
	return(obj)
}

print.ZIE <- function (x, ...)
{
	SubClass <- class(x)[1]
	cat("A", getTemp("ZIname"),
		"Import/Export definition object of subclass:", SubClass, "\n")
	cat("\n", x$description, "\n\n")
	cat("Title:  ", x$title, "\n")
	cat("Filter: ", x$filter, "\n")
	cat("Pattern:", x$pattern, "\n")
	cat("Command:", x$command, "\n")
	cat("Author: ", x$author, "\n")
	cat("Version:", x$version, "\n")
	cat("Date:    ", x$date, "\n")
	cat("License:", x$license, "\n")
	cat("Depends:", x$depends, "\n")
	cat("URL:    ", x$url, "\n")
	return(invisible(x))
}

## Import plain .tif files, with manual creation of associated .zim files
.ZIEimportTif <- ZIE(
	title       = "Tiff image files (*.tif)",
	filter      = "*.tif",
	description = c("Manual creation of ZooImage Metadata files (.zim)",
				    "given a list of directly usable TIFF images",
				    "that is, no conversion required and image names",
				    "already follow the ZooImage convention"),
	pattern     = "\\.[tT][iI][fF]$",
	command     = "zimMake(dir = Dir, pattern = Pattern, images = Files, show.log = TRUE)",
	author      = "Philippe Grosjean (phgrosjean@sciviews.org)",
	version     = "1.1-0",
	date        = "2007-02-20",
	license     = "GPL 2 or above",
	url         = "",
	depends     = "R (>= 2.4.0), zooimage (>= 1.1-0)",
	type        = "import")
 
## Import plain .jpg files, with manual creation of associated .zim files
.ZIEimportJpg <- ZIE(
	title       = "Jpeg image files (*.jpg)",
	filter      = "*.jpg",
	description = c("Manual creation of ZooImage Metadata files (.zim)",
				    "given a list of directly usable JPEG images",
				    "that is, no conversion required and image names",
				    "already follow the ZooImage convention"),
	pattern     = "\\.[jJ][pP][gG]$",
	command     = "zimMake(dir = Dir, pattern = Pattern, images = Files, show.log = TRUE)",
	author      = "Philippe Grosjean (phgrosjean@sciviews.org)",
	version     = "1.1-0",
	date        = "2007-02-20",
	license     = "GPL 2 or above",
	url         = "",
	depends     = "R (>= 2.4.0), zooimage (>= 1.1-0)",
	type        = "import")

## Complex import of images (conversion, renaming, etc.) with automatic creation
## of associated .zim files using a .zie file
.ZIEimportZie <- ZIE(
	title       = "ZooImage Import Extension (Import_*.zie)",
	filter      = "Import_*.zie",
	description = c("Run a ZIE import specification file to convert",
				    "and/or rename images and automatically create",
				    "associated ZIM files (ZooImage Metadata)"),
	pattern     = "\\.[zZ][iI][eE]$",
	command     = "zieMake(path = Dir, Filemap = Files[1], check = TRUE))",
	author      = "Philippe Grosjean (phgrosjean@sciviews.org)",
	version     = "1.1-0",
	date        = "2007-02-20",
	license     = "GPL 2 or above",
	url         = "",
	depends     = "R (>= 2.4.0), zooimage (>= 1.1-0)",
	type        = "import")

## Compile a .zie file from TemplateImport.zie and a table.txt and then compute it
.ZIEimportTable <- ZIE(
	title       = "Table and ImportTemplate.zie (*.txt)",
	filter      = "*.txt",
	description = c("Create a ZIE file by interpretting a table,",
				    "using a template file in the same directory",
				    "and named 'ImportTemplate.zie'. The resulting",
				    "ZIE file is then run to make images + metadata"),
	pattern     = "\\.[tT][xX][tT]$",
	command     = "zieCompile(path = Dir, Tablefile = Files[1], make.it = TRUE))",
	author      = "Philippe Grosjean (phgrosjean@sciviews.org)",
	version     = "1.1-0",
	date        = "2007-02-20",
	license     = "GPL 2 or above",
	url         = "",
	depends     = "R (>= 2.4.0), zooimage (>= 1.1-0)",
	type        = "import")

## Read most important EXIF data from a Digicam RAW file
.readExifRaw <- function (rawfile, full = FALSE, check = TRUE)
{	
	## Make sure dc_raw is available and rawfile exists
	if (!checkFileExists(rawfile)) return(NULL)
	
	## Temporary change directory to the one where the file is located
	filedir <- dirname(rawfile)
	if (filedir != ".") {
		inidir <- getwd()
		setwd(filedir)
		on.exit(setwd(inidir))
		rawfile <- basename(rawfile)
	}
	
	temp <- "exifdata.txt"
####	misc_dcraw(rawfile, '-i -v ', temp) 
	if (!checkFileExists(temp, message = "Error while reading exif data for '%s'"))
		return(NULL)
	
	res <- scan(temp, character(), sep = "\n", quiet = TRUE)
	if (length(res) < 6)
		return("Error getting EXIF data from '", rawfile, "'")
	
	## We replace ": " with "="
	res <- sub(": ", "=", res)
	
	## We replace all spaces by '_' (except for Filename and Timestamp,
	## first two lines!)
	res[-2] <- gsub(" ", "_", res[-2])
	if (full) {
		## Rewrite date time into yyyy-mm-dd hh:mm:ss
		datetime <- sub("^Timestamp=", "", res[2])
		lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
		newdate <- as.character(as.Date(datetime,
			format = "%a %b %d %H:%M:%S %Y"))
		Sys.setlocale("LC_TIME", lct)
		newtime <- sub("^.* (.*) [0-9]{4}$", "\\1",
			"Wed Jul 12 09:45:38 2006")
		res[2] <- paste("Timestamp=", newdate, " ", newtime, sep = "")
	} else { # Keep only most important Exif data
	    res <- res[3:7]
	}
	unlink(temp)
	return(res)
}
## example:
## setwd("g:/zooplankton/Madagascar2Macro") # Directory with the example dataset
## (Res <- .readExifRaw("Image_0742.CR2"))

## Make a comparison of two exif datasets on sensible entries
.compareExif <- function (Exif1, Exif2)
{
	dif <- character(0)
	## Need same 'Camera', 'ISO_speed', 'Shutter', 'Aperture', 'Focal_Length'
	### TODO: make it work for larger Exif dataset. Currently requires that the 
	###       fields are restricted to strict equal data
	if (length(Exif1) != length(Exif2)) {
	    dif <- "Not same size for both Exif data!"
	} else {
	    difpos <- sort(Exif1) != sort(Exif2)
	    if (any(difpos)) dif <- "Exif data are not identical!"
	}
	return(dif)
}

.isTestFile <- function (File)
{
	## Determine if a given file is a test file (a file with first line being
	## 'ZI1est' and with size < 1000)
	if (file.info(File)$size > 1000) return(FALSE)
	checkFirstLine(File, "ZItest")
}

## Check a blank-field image, either in .pgm or in .tif format	
.checkBF <- function (BFfile)
{	
	if (!checkFileExists(BFfile, message = "Blank-field file '%s' not found!"))
		return(NULL)

	## Is it a test file?
	if (.isTestFile(BFfile))
		return(character(0)) # We behave like if the file was correct!
	
	msg <- character(0)
	filedir <- dirname(BFfile)
	if (filedir != ".") {
		## Temporary change directory to the one where the file is located
		inidir <- getwd()
		setwd(filedir)
		on.exit(setwd(inidir))
		BFfile <- basename(BFfile)
	}
	
	## The command to use depends on the format of the image (determined on the
	## extension)
	ext <- tolower(rev(strsplit(BFfile, "\\.")[[1]])[1])
	pgmfile <- BFfile
	if (ext == "tif") {
		## First, convert into a .pgm file
		pgmfile <- paste(BFfile, "pgm", sep = ".")
####		netpbm_tifftopnm( BFfile, pgmfile )
		delfile <- TRUE
		ext <- "pgm"
	} else delfile <- FALSE
	
	if (ext != "pgm")
		return(paste("Unrecognized image format for '", BFfile, "'", sep = ""))
	
####	BF <- netpbm_pgmhist(pgmfile, delete = delfile)
	
	## Make sure we work with 16bit images
	if (max(BF$Gray) < 256) {
		msg <- c(msg, "Blank-field seems to be a 8bit image (16bit required)")	
	} else {
		## Look at darkest value with at least 10 points
		BF <- BF[BF$Count >= 10, ]
		darkpart <- min(BF$Count)
		
		## Eliminate values with low number of points
		BF <- BF[BF$Count >= 100, ]
		
		## Check range for these values
		rngBF <- range(BF$Gray)
		if (rngBF[2] > 65500)
			msg <- c(msg, "Blank-field is overexposed")
		if (rngBF[2] < 60000)
			msg <- c(msg, "Blank-field is underexposed or contains too dark areas")
		if ((rngBF[2] - rngBF[1]) > 15000)
			msg <- c(msg, "Blank-field is too heterogeneous")
		if ((rngBF[1] - darkpart) > 15000)
			msg <- c(msg, "Blank-field contains dark zones (dust?)")
	}
	return(msg)
}
## example:
## setwd("g:/zooplankton/Madagascar2Macro") # Directory with the example dataset
## .checkBF("test.pgm")
## .checkBF("test.tif")

## Make a blank-field correction on File, given a BFfile (blank-field) 
## Both files must be 16bit gray PGM images
## The resulting file has same name as File, but with a .tif extension instead
## of .pgm
## The function returns TRUE in case of success... or an explicit error message
.BFcorrection <- function (File, BFfile, deleteBF = TRUE, check = TRUE)
{
	on.exit({
		unlink(imgFile)
		if (deleteBF) unlink(imgBFfile)
	})
	if (!checkFileExists(File, "pgm")) return(NULL)
	if (!checkFileExists(BFfile, "pgm", message = "Blank-field file '%s' not found"))
		return(NULL)
	
	## Check that the various scripts are available
	#checkCapable("pnm2biff")
	#checkCapable("statistics")
	#checkCapable("divide")
	#checkCapable("biff2tiff")
	
	## Switch to the directory of File
	filedir <- dirname(File)
	if (filedir != ".") {
		## Temporary change directory to the one where the file is located
		inidir <- getwd()
		setwd(filedir)
		on.exit(setwd(inidir), add = TRUE)
		File <- basename(File)
	}
	
	## Determine the name of the various files
	fileNoExt <- noExtension(File)
	imgFile <- paste(fileNoExt, "img", sep = ".")
	imgcorrFile <- paste(fileNoExt, "coor.img", sep = "")
	tifFile <- paste(fileNoExt, "tif", sep = ".")
	imgBFfile <- paste(noExtension(BFfile), "img", sep = ".")

    ## Is File a test file?
	if (.isTestFile(File)) {
		## We behave like if the file was corrected, but just copy the content
		## of File into tifFile
		file.copy(File, tifFile)
		
		## Simulate creation of the .img blank-field
		if (!deleteBF) file.copy(BFfile, imgBFfile)
		return(TRUE)
	}

	## Convert PGM files into BIFF
####	xite_pnm2biff(File, imgFile)
####	xite_pnm2biff(BFfile, imgBFfile)
	
	## Get the mean gray level of the blank-field
####	meangray <- xite_statistics(imgBFfile)
	
	## Eliminate the blank field
####	res <- xite_divide(meangray, imgFile, imgBFfile, imgcorrFile)
	
	## Make the tiff file
####	xite_biff2tiff(imgcorrFile, tifFile)
	
	return(TRUE) # Everything is fine!
}
## example: 
## setwd("g:/zooplankton/madagascar2macro")
## .BFcorrection("_CalibOD03.pgm", "_CalibBF03.pgm")

## Convert a RAW file (digital camera) into a pgm file
### TODO: can we not do this in the the ImageJ plugin directly
.rawConvert <- function (RawFile, OutputFile = "fileconv.pgm",
DcRawArgs = "-v -c -4 -q 3 -t 0 -k 0", fake = FALSE, replace = FALSE,
check = TRUE)
{	
	## Check if the output file already exists
	if (file.exists(OutputFile)) {
		## If we want to replace existing file, delete it, otherwise, we are done
		if (replace) unlink(OutputFile) else return(TRUE)
	}
	
	## Check if RawFile exists
	if (!checkFileExists(RawFile)) return(FALSE)
	
	## Do a fake convert
	if (fake) { # Create a test file with just ZItest in it
		cat("ZItest\n", file = OutputFile)
		return(TRUE)
	}
	
	## Do the conversion using dc_raw
	## Check that the system is capable of doing the conversion
	if (check) {
		#checkCapable("dc_raw")
		#checkCapable("ppmtopgm")
	}
	
	if (isWin()) {
		## Convert the RAW file into PPM file (48bit color)
		## We have to do it in two steps because windows lack of proper piping
####		misc_dcraw(RawFile, DcRawArgs , "RAWTEMP.PPM")

		## Convert from 48bit color to 16bit grayscale
####		netpbm_ppmtopgm("RAWTEMP.PPM", OutputFile)
	} else {
		## One step conversion (no tempfile)
		cmd <- sprintf('dcraw %s %s | ppmtopgm > "%s"' , 
			DcRawArgs, RawFile, OutputFile)
		res <- try(system(cmd), silent = TRUE)
		if (!checkFileExists(OutputFile, message = "Error while converting"))
			return(FALSE)
	}
	
	## Everything was fine
	return(TRUE)
}
## example:
## setwd("d:/ZI examples/MacroG16-example")
## .rawConvert("Image_3822.CR2", fake = TRUE)
## .rawConvert("Image_3822.CR2")
