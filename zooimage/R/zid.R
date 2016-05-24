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

## Check consistency of a zooimage directory before creating .zid or .zidb file
zidVerify <- function (zidir, type = c("ZI1", "ZI2", "ZI3"),
check.vignettes = TRUE)
{	
	## Check the format of the file
	## This should be a directory containing XXX+YY_dat1.zim files
	## + .jpg or .png files (vignettes)
	if (any(!type %in% c("ZI1", "ZI2", "ZI3"))) {
		warning("only 'ZI1', 'ZI2' or 'ZI3' are currently supported for 'type'")
		return(invisible(FALSE))
	}
	
	## Check the list of _dat1.zim
	dat1files <- zimDatList(zidir)
	if (!length(dat1files)) {
		warning("no '_dat1.zim' file!")
		return(invisible(FALSE))
	}
	
    ## Check the content of all these "_dat1.zim" files 
	## and retrieve the number of items measured
	dat1files <- sort(dat1files)
	## Default to -1 for corrupted dat1 files
	nitems <- sapply(dat1files, function(x) {
		zimVerify(file.path(zidir, x))
	})
	ok <- all(nitems >= 0)
	if (!ok) {
		warning("corrupted '_dat1.zim' files: ", paste(dat1files[nitems < -1],
			collapse = ", "))
		return(invisible(FALSE))
	}
	
	## Check the vignettes
	if (isTRUE(as.logical(check.vignettes))) {
        ## Check that we have corresponding vignettes (XXX+YY_ZZZ.jpg/png files)
    	samples <- sub("_dat1[.]zim$", "", dat1files)
		
    	## Check the content of the directory for .jpg or .png files
    	for (i in 1:length(samples)) {		
			## List the jpegs
    		regex <- gsub("[+]", "[+]", samples[i])
    		regex <- gsub("[.]", "[.]", regex)
    		regex2 <-  paste("^", regex, "_[0-9]+[.]jpg$", sep = "")
			vigstype <- "jpg"
    		vigs <- dir(zidir, pattern = regex2)
			if (!length(vigs)) { # Try also for .png vignettes
				regex2 <-  paste("^", regex, "_[0-9]+[.]png$", sep = "")
				vigstype <- "png"
				vigs <- dir(zidir, pattern = regex2)
			}
			
    		## Get their numbers, sort them, and make sure none is missing
    		n <- nitems[i]
			## If impossible to know how many items, just count vignettes
    		if (n < 1) n <- length(vigs)
    		
			## Construct a vector with names of vignettes as they should be
    		chkvigs <- paste(samples[i], "_", 1:n, ".", vigstype, sep = "")
    		if (length(vigs) == 0 && length(chkvigs) > 0) {
				warning("no vignettes for ", samples[i])
				ok <- FALSE
            } else if (length(chkvigs) != length(vigs) ||
				!all(sort(chkvigs) == sort(vigs))) {
				warning("mismatch vignettes for ", samples[i])
				ok <- FALSE 
			}
        }
    } 

	invisible(ok)
}

zidVerifyAll <- function (path = ".", samples = NULL,
type = c("ZI1", "ZI2", "ZI3"), check.vignettes = TRUE)
{	
	## Verify all of these directories
	if (any(!type %in% c("ZI1", "ZI2", "ZI3"))) {
		warning("only 'ZI1', 'ZI2' or 'ZI3' are currently supported for 'type'")
		return(invisible(FALSE))
	}
	
	## First, switch to that directory
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))	
	path <- "."	# Indicate we are now in the right path
	
	## Process the list of samples
	if (!length(samples)) {	# Compute them from path
		d <- dir(path, pattern = "^[^_]")	# All items not starting with '_'
		samples <- unique(d[file.info(d)$isdir])	# Keep only directories
	}
	
	## If there is no dir, exit now
	if (!length(samples)) {
		warning("There is no directories to verify in ", getwd())
		return(invisible(FALSE))	
	}
	
	## Start the process
	smax <- length(samples)
	message("Verification of .zid content...")
	flush.console()
	ok <- batch(samples, zidVerify, type = type,
		check.vignettes = check.vignettes, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(samples),
			" samples pass verification (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Compress one sample as a single .zid zipped file	
zidCompress <- function (zidir, type = c("ZI1", "ZI2", "ZI3"), check = TRUE,
check.vignettes = TRUE, replace = FALSE, delete.source = replace)
{		
	## Check the format
	if (any(!type %in% c("ZI1", "ZI2", "ZI3"))) {
		warning("only 'ZI1', 'ZI2' or 'ZI3' are currently supported for 'type'")
		return(invisible(FALSE))	
	}
	
	## We need to switch to the root of sample dir first for correct path
	## in the zip file
	rootdir <- dirname(zidir)
	initdir <- setwd(rootdir)
	on.exit(setwd(initdir))
	zidir <- basename(zidir) # Use only the latest dir (the "sample dir")
	
	## The .zid file is located in the "root" dir, same name as the
	## "sample dir", with .zid extension
	zidfile <- paste(zidir, "zid", sep = ".")
	if (!isTRUE(as.logical(replace)) && file.exists(zidfile)) {
		## It is not advised to delete source without rebuilding the .zid file
		## but it was expressly asked!
		### TODO: verify we have the same files in the .zid and initial dir
		## before deleting files!
		if (delete.source && file.exists(zidir))
			unlink(zidir, recursive = TRUE)
		return(invisible(TRUE))	# Nothing else to do
	}
	
	## Make sure everything is fine for this directory
	if (isTRUE(as.logical(check)))
		if (!zidVerify(zidir, type = type, check.vignettes = check.vignettes))
			return(invisible(FALSE))
	
	## Make sure the .RData file is created (or refreshed)
	if (!zidDatMake(zidir, type = type, replace = replace))
		return(NULL)
	
	## Do compress the directory in the .zip file
	## Copy or move all corresponding files to a .zid zip-compressed file
	res <- zip(zidfile, zidir, flags = "-rq9X")
	
	## Do we delete sources?
	if (isTRUE(as.logical(delete.source)))
		unlink(zidir, recursive = TRUE)
	
	invisible(res != 0)
}

## Compress all data in the corresponding directory
zidCompressAll <- function (path = ".", samples = NULL,
type = c("ZI1", "ZI2", "ZI3"), check = TRUE, check.vignettes = TRUE,
replace = FALSE, delete.source = replace)
{ 
	if (any(!type %in% c("ZI1", "ZI2", "ZI3"))) {
		warning("only 'ZI1', 'ZI2' or 'ZI3' are currently supported for 'type'")
		return(invisible(FALSE))
	}
	
	## First, switch to that directory                                       
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- "."	# Indicate we are now in the right path
	
	## Get the list of samples to process
	if (!length(samples)) {	# Compute them from path
		d <- dir(path, pattern = "^[^_]")	# All items not starting with '_'
		samples <- unique(d[file.info(d)$isdir])	# Keep only directories
	}
	
	## If there is no dir, exit now
	if (!length(samples)) {
		warning("there is no directories to process in ", getwd())
		return(invisible(FALSE))	
	}
		
	## Start the process
	if (isTRUE(as.logical(check)))
		if (!zidVerifyAll(path = path, samples = samples, 
			check.vignettes = check.vignettes))
			return(invisible(FALSE))
		
	## Compress these files
	message("Compression of .zid data...")
	flush.console()
	ok <- batch(samples, zidCompress, type = type, check = FALSE, 
			check.vignettes = check.vignettes, replace = replace,
			delete.source = delete.source, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(samples),
			" items were compressed in .zid files (see .last.batch)")
		invisible(FALSE)
	} else {
		## Possibly clean the whole directory (move .zim files to \_raw
		## and delete the \_work subdir if everything is fine
		zidClean(path = path, samples = samples)	
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Clean Zid (eliminate the _work subdirectory and move initial data to _raw)
zidClean <- function (path = ".", samples = NULL)
{
	## Do we have samples to process
    if (!length(samples)) return(invisible(FALSE))
	
    ## First, switch to that directory
    if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	    
	## Identify paths
	message("Cleaning directory...")
	zimfiles   <- zimList( "." )
	zimsamples <- sub("^(.*)[+].+", "\\1", zimfiles)
	
    ## Keep only those .zim files related to samples
    zimfiles <- zimfiles[zimsamples %in% samples]
	
	## Process
    if (length(zimfiles)) {
        rawdir <- file.path(".", "_raw")
        
		## If the _raw subdirectory does not exists, create it
        if (!file.exists(rawdir)) dir.create(rawdir)
		copyto <- file.path(".", "_raw", zimfiles)
        
		## Move these .zim files
        for (i in 1:length(zimfiles))
			file.rename(zimfiles[i], copyto[i])
    }
	
    ## Delete completely the _work subdirectory
    unlink(file.path(".", "_work"), recursive = TRUE)
	
	return(invisible(TRUE))
}

## Uncompress a .zid file to get all its content. 
## Use 'delete.source = TRUE' with caution!
zidUncompress <- function (zidfile, path = dirname(zidfile),
delete.source = FALSE)
{	
	## Check if the file provided is a .zid file, and if it exists
	if (!checkFileExists(zidfile, extension = "zid"))
		return(invisible(FALSE))
	
	message("Unzipping '", zidfile, "' ...")
	
	## Uncompress it
	if (!length(tryCatch(unzip(zidfile, overwrite = FALSE, exdir = "."),
			error = function (e) warning(e),
			warning = function (w) return()))) {
		message("    ... not done!")
		return(invisible(FALSE))
	}
	
	## Do we delete sources?
	if (isTRUE(as.logical(delete.source))) unlink(zidfile)
	
	## Invisibly indicate success
	invisible(TRUE)
}

## Uncompress all .zid files in the 'path.extract' directory
zidUncompressAll <- function (path = ".", zidfiles = zidList(path),
path.extract = path, skip.existing.dirs = TRUE, delete.source = FALSE)
{	
	## Initial checks
	if (!length(zidfiles)) {
        warning("no ZID files!")
		return(invisible(FALSE))
	}
	
	## Start the process
	ok <- TRUE

	## Check that dirs / files with corresponding names exist in path.extract
	checkdirs  <- file.path(path.extract, noExtension(zidfiles))
	fileExists <- file.exists(checkdirs) & !file.info(checkdirs)$isdir
	dirExists  <- file.exists(checkdirs) & file.info(checkdirs)$isdir
	
	## If any file not being a dir exist there, stop the process
	if (any(fileExists)) {
        warning("one or several files have same name as uncompressed dirs!")
		return(invisible(FALSE))
	}
	
	## Should we eliminate files whose corresponding dirs exist?
	if (skip.existing.dirs && any(dirExists)) {
		cat(sum(dirExists), "file(s) already uncompressed skipped!\n")
        warning(paste("Skipping already uncompressed file(s):", 
			paste(zidfiles[dirExists], collapse = ",")))
	}
	zidfiles <- zidfiles[!dirExists]
	
	## Decompress the files remaining in the list
	smax <- length(zidfiles)
	if (!length(zidfiles)) {
		message("-- Done! - (nothing to decompress)")
		return(invisible(TRUE))
	}
	
	## Uncompress these files
	message("Decompression of ZID archives...")
	flush.console()
	ok <- batch(zidfiles, zidUncompress, path = path.extract,
		delete.source = delete.source, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(zidfiles),
			" ZID files were uncompressed (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}	
}

## Make a .RData file that collates together data from all the "_dat1.zim" files
## of a given sample
zidDatMake <- function (zidir, type = "ZI3", replace = FALSE) 
{
    if (any(!type %in% c("ZI1", "ZI2", "ZI3"))) {
		warning("only 'ZI1', 'ZI2' or 'ZI3' are currently supported for 'type'")
		return(invisible(FALSE))
	}
		
    RDataFile <- file.path(zidir, paste0(basename(zidir), "_dat1.RData"))
    
    ## File already exists
    if (file.exists(RDataFile) && !replace) 
        return(invisible(TRUE))
    
	ok <- TRUE
    dat1files <- zimDatList(zidir)

    ## Create _dat1.zim file if it is missing (for FlowCAM data)
    if (!length(dat1files)) {
        SmpDir <- dirname(zidir) 
        zimDatMakeFlowCAM(file.path(SmpDir,
			paste(basename(zidir), "zim", sep = ".")))
        dat1files <- zimDatList(zidir)
        if (!length(dat1files)) {
            warning("no '_dat1.zim' file!")
			return(invisible(FALSE))
		}
    }
    
    dat1files <- sort(dat1files)
    fractions <- sampleInfo(dat1files, "fraction")

    ## Avoid collecting duplicate informations about fractions
    fracdup <- duplicated(fractions)
    results <- lapply(seq.int(1, length(dat1files)), function (i) {
        dat1path <- file.path(zidir, dat1files[i])
        if (!isZim(dat1path)) return(invisible(FALSE))
        
        ## Read the header
        Lines <- scan(dat1path, character(), sep = "\t", skip = 1, 
            blank.lines.skip = FALSE, flush = TRUE, quiet = TRUE, 
            comment.char = "#")
        if (length(Lines) < 1) {
            warning( dat1files[i], " is empty, or is corrupted")
            return(invisible(FALSE))
        }
        
        ## Trim leading and trailing spaces in Lines
        Lines <- trimString(Lines)
        
        ## Convert underscore to space
        Lines <- underscoreToSpace(Lines)
        
        ## Determine start of the measurements table (it is '[Data]' header)
        endhead <- tail(which(Lines == "[Data]"), 1)
        if (!is.null(endhead) && endhead > 1) 
            Lines <- Lines[seq.int(1, endhead - 1)]
        
        ## Decrypt all lines, that is, split on first occurrence
		## of "=" into 'tag', 'value' and separate into sections
        if (!is.null(Lines))
            meta <- parseIni(Lines, sub("_dat1[.]zim$", "", fractions[i]))

        if (!is.null(endhead)) {
            mes <- read.table(dat1path, header = TRUE, sep = "\t", 
                dec = ".", as.is = FALSE, skip = endhead + 1, 
                comment.char = "#", na.strings = "null")
            ## We have several problems here:
            ## 1) There is sometimes a column full of NAs at the end.
            ##    This is because ImageJ adds an extra tab at the end of the line.
            ## [RF] FIXME: this should not be the case anymore because we have
            ## more control on what ImageJ is doing
			## [PhG] We keep this here anyway for old datasets!
            if (all(is.na(mes[, ncol(mes)]))) 
                mes <- mes[, -ncol(mes)]
            ## 2) The first column is the 'Item', but its name '!Item' is
            ##    transformed into 'X.Item'
            ## 3) The '%Area' is transformed into 'X.Area'
            Names <- names(mes)
            if (Names[1] == "X.Item") Names[1] <- "Item"
            if ("X.Area" %in% Names) Names[Names == "X.Area"] <- "PArea"
            ## Invert 'Item' and 'Label'
            mes <- mes[, c(2, 1, 3:ncol(mes))]
            Names <- Names[c(2, 1, 3:length(Names))]
            names(mes) <- make.names(Names, unique = TRUE)
            Sub <- meta$Subsample
            Sub$Dil <- 1/(Sub$SubPart * Sub$CellPart * Sub$Replicates * 
                Sub$VolIni)
            mes$Dil <- rep(Sub$Dil[Sub$Label == fractions[i]], nrow(mes))
        } else {
            mes <- NULL
        }
        list(meta = meta, mes = mes)
    })
    notnull.filter <- Negate(is.null)
    results <- Filter(notnull.filter, results)
    list.allmeta <- Filter(notnull.filter, lapply(results, "[[", "meta"))
    list.allmes <- Filter(notnull.filter, lapply(results, "[[", "mes"))
    
    listCombine <- function (.list) {
        force(.list)
        mergefun <- function (x, y) {
            if (all(sort(names(x)) == sort(names(y)))) {
                rbind(x, y)
            } else {
                merge(x, y, all = TRUE)
            }
        }
        Reduce(mergefun, .list)
    }
	
	listMerge <- function (x, y) {
		xitems <- names(x)
		yitems <- names(y)
		xandy <- xitems[xitems %in% yitems]
		xonly <- xitems[!(xitems %in% xandy)]
		yonly <- yitems[!(yitems %in% xandy)]
		## First merge common items
		if (length(xandy) > 0) {
			res <- lapply(xandy, function (item) {
				merge(x[[item]], y[[item]], all = TRUE)
			})
			names(res) <- xandy
		} else {
			res <- list()
		}
		if (length(xonly) > 0) res[xonly] <- x[xonly]
		if (length(yonly) > 0) res[yonly] <- y[yonly]
		res
	}

    list.allmeta <- list.allmeta[!fracdup] # only the levels of not duplicated metadata
    lmeta <- length(list.allmeta[])
    if (lmeta == 1) {
        allmeta <- listCombine(list.allmeta)
    } else {
        allmeta <- NULL
        for (i in 1:(lmeta - 1))
            allmeta <- listMerge(list.allmeta[[i]], list.allmeta[[i + 1]])
    }

    allmes <- listCombine(list.allmes)
    rownames(allmes) <- 1:nrow(allmes)
    Names <- names(allmes)
    
	## Calculate an ECD from Area if there is not one yet
    if (!"ECD" %in% Names && "Area" %in% Names) {
        ECD <- ecd(allmes$Area)
        allmes <- data.frame(allmes[, 1:2], ECD = ECD, allmes[, 3:ncol(allmes)])
    }
    attr(allmes, "metadata") <- allmeta
    class(allmes) <- c("ZI3Dat", "ZIDat", "data.frame")
    ZI.sample <- allmes
    save(ZI.sample, file = RDataFile, ascii = FALSE, version = 2,
		compress = TRUE)
    if (ok) ok <- file.exists(RDataFile)
	if (!ok) warning("problem creating the RData file")
    return(invisible(ok))
}

## Read the .Rdata in a .zid file or corresponding directory
zidDatRead <- function (zidfile)
{	
	## Identify the file and stop if it does not exists
	sample <- noExtension(zidfile)
	RdataFile <- paste0(sample, "_dat1.RData")
	deletefile <- FALSE
	if (!checkFileExists(zidfile, message = "%s not found!")) return(NULL)
	
	## Treat different kind of files
	if (!hasExtension(zidfile, "zid")) {
		# Is it a directory?
		if (file.info(zidfile)$isdir) {
			# Is there a .RData file in this directory?
			rdata <- file.path(zidfile, RdataFile)
			if (!file.exists(rdata)) {
				# Try to create it
				zidDatMake(zidfile)
				if (!checkFileExists(rdata,
					message = "error creating the RData file"))
					return(NULL)
			}
		} else {
			warning("unrecognized file: ", zidfile)
			return(NULL)
		}
	} else {  # This is a .zid file
		rdata <- file.path(sample, RdataFile)
		
		
		zidExtract <- function (file, zidfile) {
			tmpd <- tempdir()
			unzip(zidfile, file, exdir = tmpd, overwrite = TRUE,
				junkpaths = TRUE)
			res <- file.path(tmpd, file)
			if (file.exists(res)) res else NULL
		}		
		rdata <- zidExtract(rdata, zidfile)
		if (!length(rdata)) {
			warning("error reading RData file from ", basename(zidfile))
			return(NULL)
		}
		deletefile <- TRUE
	}
	
	## Load that file
	ZI.sample <- NULL
	load(rdata)
	
	## Delete the file
	if (deletefile) {
		unlink(rdata)
		# If the directory is empty, delete it also
		datadir <- file.path(tempdir(), sample)
		if (file.exists(datadir) && !length(dir(datadir)))
			unlink(datadir)
	}
	
	## Set the class 
	if (!inherits(ZI.sample, "ZIDat") && inherits(ZI.sample, "data.frame"))
		class(ZI.sample) <- c("ZI3Dat", "ZIDat", "data.frame")
	return(ZI.sample)
}
