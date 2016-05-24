## Copyright (c) 2012, Ph. Grosjean <phgrosjean@sciviews.org> & K. Denis
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

## Make a ZooImage database file for one sample
zidbMake <- function (zidir, zidbfile = paste0(zidir, ".zidb"),
zisfile = file.path(dirname(zidir), "Description.zis"), type = "ZI3",
check = TRUE, check.vignettes = TRUE, replace = FALSE,
delete.source = replace)
{		
	## Check the format
	if (type != "ZI3") {
		warning("only 'ZI3' is currently supported for 'type'")
		return(invisible(FALSE))	
	}
	
	if (!isTRUE(as.logical(replace)) && file.exists(zidbfile)) {
		## Nothing to do... the file already exists
		if (isTRUE(as.logical(delete.source)) &&
			file.exists(zidir) && file.info(zidir)$isdir)
			unlink(zidir, recursive = TRUE)
		return(invisible(TRUE))	# Nothing else to do
	}
	
	## Make sure everything is fine for this directory
	if (isTRUE(as.logical(check)))
		if (!zidVerify(zidir, type = type, check.vignettes = check.vignettes))
			return(invisible(FALSE))
	
	## Make sure the .RData file is created (or refreshed)
	if (!zidDatMake(zidir, type = type, replace = replace))
		return(invisible(FALSE))
	
    ## List all vignettes
    Vigs <- dir(zidir, pattern = "\\.jpg$", full.names = TRUE)
	## Maybe the vignettes are in .png format...
	if (!length(Vigs)) {
		Vigs <- dir(zidir, pattern = "\\.png$", full.names = TRUE)
		VigType <- "png"
	} else VigType <- "jpg"
	if (!length(Vigs)) {
		warning("No vignettes found (JPEG or PNG files)")
		return(invisible(FALSE))
	}
    ## List all .zim files
    Zims <- dir(zidir, pattern = "\\.zim$", full.names = TRUE)
	if (!length(Zims)) {
		warning("No ZIM files found!")
		return(invisible(FALSE))
	}
    
	## Make sure data from the .zis file are correct
	if (!checkFileExists(zisfile, "zis", force.file = TRUE))
		return(invisible(FALSE))
	zisData <- zisRead(zisfile)
	isSample <- (zisData$Label == basename(zidir))
	if (!length(isSample) || sum(isSample) < 1) {
		warning("Incorrect .zs file, or the file does not contain data for this sample")
		return(invisible(FALSE))
	}
    ## Extract data for this sample
	zisData <- zisData[isSample, ]
	## TODO: may be check that a minimum set of variables is there...
	
    ## Create the .zidb file: put all vignettes there, plus the .RData file
	message("Creating the ZIDB file...")
    filehashOption(defaultType = "DB1")
    unlink(zidbfile)
    dbCreate(zidbfile)
    db <- dbInit(zidbfile)

    ## Indicate which zooimage version and which image type we use
    dbInsert(db, ".ZI", 3)
    dbInsert(db, ".ImageType", VigType)

    ## Read each vignette in turn and add it to the database
	message("Adding vignettes to ZIDB file...")
    VigExt <- paste0("\\.", VigType, "$")
	for (i in 1:length(Vigs)) {
    	Vig <- Vigs[i]
    	VigName <- sub(VigExt, "", basename(Vig))
    	VigSize <- file.info(Vig)$size
		if (is.na(VigSize)) {
			warning("file '", Vig, "' not found, or of null length")
			return(invisible(FALSE))
		}
    	dbInsert(db, VigName, readBin(Vig, "raw", VigSize + 100))
    }
    
    ## Add .zim files to db
	message("Adding data from ZIM files to ZIDB file...")
    for (i in 1:length(Zims)) {
    	Zim <- Zims[i]
    	ZimName <- sub("\\.zim$", "", basename(Zim))
    	ZimSize <- file.info(Zim)$size
		if (is.na(ZimSize)) {
			warning("file '", Zim, "' not found or of null length")
			return(invisible(FALSE))	
		}
    	dbInsert(db, ZimName, readBin(Zim, "raw", ZimSize + 100))
    }

    ## Add zis info to db
	message("Adding sample data to ZIDB file...")
	dbInsert(db, ".SampleData", zisData)

    ## Add the data frame with all data and metadata to the file
	message("Adding R data to ZIDB file...")
    zidat <- file.path(zidir, paste0(basename(zidir), "_dat1.RData"))
    obj <- load(zidat)
	if (length(obj) != 1) {
		warning("Error loading ", zidat)
		return(invisible(FALSE))
	}
    dbInsert(db, ".Data", get(obj))

	## Do we delete sources?
    if (isTRUE(as.logical(delete.source)))
        unlink(zidir, recursive = TRUE)
		
	message("-- Done! --")
	
	## Indicate success...
	invisible(TRUE)
}

## Make all .zidb files for data in the corresponding directory
zidbMakeAll <- function (path = ".", samples,
zisfiles = file.path(dirname(samples), "Description.zis"), type = "ZI3",
check = TRUE, check.vignettes = TRUE, replace = FALSE, delete.source = replace)
{ 
	if (type != "ZI3")
		stop("only 'ZI3' is currently supported for 'type'")
	
	## First, switch to that directory                                       
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- "."	# Indicate we are now in the right path
	
	## Get the list of samples to process
	if (missing(samples) || !length(samples)) {	# Compute them from path
		## All dirs not starting with '_'
		dirs <- dir(path, pattern = "^[^_]", full.names = TRUE)
		samples <- unique(dirs[file.info(dirs)$isdir]) # Keep only directories
	}
	
	## If there is no dir, exit now
	if (!length(samples)) {
		warning("there are no directories to process in ", getwd())
		return(invisible(FALSE))
	}
	
	## Check zisfiles and make sure the vector has same length as samples
	## possibly recycling the file(s)
	zisfiles <- as.character(zisfiles)
	if (!length(zisfiles)) {
		warning("You must provide at least one ZIS file with samples characteristics")
		return(invisible(FALSE))
	}
	if (!checkFileExists(zisfiles, "zis", force.file = TRUE))
		return(invisible(FALSE))
	zisfiles <- rep(zisfiles, length.out = length(samples))
		
	## Possibly verify the files
	if (isTRUE(as.logical(check)))
		if (!zidVerifyAll(path = path, samples = samples, 
			check.vignettes = check.vignettes))
			return(invisible(FALSE))
		
	## Create the .zidb files
	message("Creation of ZIDB files...")
	flush.console()
	zidbMakeOne <- function (item, samples, zisfiles, type, check.vignettes,
		replace, delete.source)
		zidbMake(samples[item], zisfile = zisfiles[item], type = type,
			check = FALSE, check.vignettes = check.vignettes, replace = replace,
			delete.source = delete.source)
	items <- 1:length(samples)
	ok <- batch(items, zidbMakeOne, samples = samples, zisfiles = zisfiles,
			type = type, check.vignettes = check.vignettes, replace = replace,
			delete.source = delete.source, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(samples),
			" items were correctly processed (see .last.batch)")
		invisible(FALSE)
	} else {
		## Possibly clean the whole directory (move .zim files to \_raw
		## and delete the \_work subdir if everything is fine
		zidClean(path = path, samples = samples)	
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Convert .zid file to .zidb file
zidToZidb <- function (zidfile, zisfile = file.path(dirname(zidfile),
"Description.zis"), replace = FALSE, delete.source = replace)
{
    if (!file.exists(paste0(zidfile, "b")) || isTRUE(as.logical(replace))) {
		ZidDir <- sub("\\.zid$", "", zidfile)
		IniDir <- dirname(zidfile)
    
		## Unzip the file...
		message("Unzipping ZID file '", basename(zidfile), "' ...")    
		if (!length(tryCatch(unzip(zidfile, overwrite = replace,
			junkpaths = FALSE, exdir = IniDir), error = function (e) warning(e),
			warning = function (w) return()))) {
			message("    ... not done!")
			return(invisible(FALSE))
		}

		## Make sure ZidDir is created...
		if (!checkDirExists(ZidDir,
			message = 'expected unzipped dir "%s" not found'))
			return(invisible(FALSE))

		## Create the .zidb file
		res <- zidbMake(zidir = ZidDir, type = "ZI3", check = TRUE,
			check.vignettes = TRUE, replace = replace,
			delete.source = delete.source)
	
	} else res <- TRUE
	
	# Do we have to delete the zidfile?
	if (res && isTRUE(as.logical(delete.source))) unlink(zidfile)
		
	message("-- Done! --")
		
	invisible(res)
}

## Convert all .zid files to .zidb files
zidToZidbAll <- function (path = ".", zidfiles, zisfiles =
file.path(dirname(zidfiles), "Description.zis"), replace = FALSE,
delete.source = replace)
{
    ## First, switch to that directory
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- "."	# Indicate we are now in the right path
	
	## Get the list of zidfiles to process
	if (missing(zidfiles) || !length(zidfiles))	# Compute them from path
		zidfiles <- dir(path, pattern = extensionPattern("zid"),
			full.names = TRUE) # All .zid files

	## If there is no zidfiles in this dir, exit now
	if (!length(zidfiles)) {
		warning("There is no ZID files to process in ", getwd())
		return(invisible(FALSE))	
	}

	## Make sure there is no path associated
	#if (!all(zidfiles == basename(zidfiles))) {
	#	warning("You cannot provide paths for ZID files, just file names")
	#	return(invisible(FALSE))
	#}
	
	## Check zisfiles and make sure the vector has same length as zidfiles
	## possibly recycling the file(s)
	zisfiles <- as.character(zisfiles)
	if (!length(zisfiles)) {
		warning("You must provide at least one ZIS file with samples characteristics")
		return(invisible(FALSE))
	}
	if (!checkFileExists(zisfiles, "zis", force.file = TRUE))
		return(invisible(FALSE))
	zisfiles <- rep(zisfiles, length.out = length(zidfiles))
	
	## Create the .zidb files from the .zid files
	message("Conversion of ZID to ZIDB files...")
	flush.console()
	zidConvertOne <- function (item, zidfiles, zisfiles, replace, delete.source)
		zidToZidb(zidfiles[item], zisfile = zisfiles[item], replace = replace,
			delete.source = delete.source)
	items <- 1:length(zidfiles)
	ok <- batch(items, zidConvertOne, zidfiles = zidfiles, zisfiles = zisfiles,
		replace = replace, delete.source = delete.source, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(zidfiles),
			" items were correctly converted (see .last.batch)")
		invisible(FALSE)
	} else {	
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Convert a .zidb file to a .zid file
zidbToZid <- function (zidbfile, zisfile = file.path(dirname(zidbfile),
"Description.zis"), replace = FALSE, delete.source = replace)
{
	zidfile <- paste(basename(zidbfile), "zid", sep = ".")
	if (!isTRUE(as.logical(replace)) && file.exists(zidfile)) {
		## It is not advised to delete source without rebuilding the .zid file
		## but it was expressly asked!
		### TODO: verify we have the same data in the .zid and .zidb files
		### before deleting the .zidb file!
		if (delete.source && file.exists(zidbfile))
			unlink(zidbfile)
		return(invisible(TRUE))	# Nothing else to do
	}
	
	if (!file.exists(zidfile) || isTRUE(as.logical(replace))) {
		ZidDir <- sub("\\.zidb$", "", zidbfile)
		## Create the directory to extract data
		dir.create(ZidDir)
		## Link database to objects in memory
		Zidb <- zidbLink(zidbfile)
		## All files in Zidb
		AllFiles <- ls(Zidb) # List vars not starting with . => zims + vignettes

		# .zim files
		isZimFile <- grep("_dat.$", AllFiles)
		ZimNames <- AllFiles[isZimFile]
		message("Extracting data from ZIM files...")
		for (ZimName in ZimNames)
		    writeBin(Zidb[[ZimName]],
				con = file.path(ZidDir, paste0(ZimName, ".zim")))
    
		## Vignettes
		VignNames <- AllFiles[-isZimFile]
		message("Extracting vignettes...")
		extension <- Zidb$.ImageType
		for(i in 1:length(VignNames)){
		    writeBin(Zidb[[VignNames[i]]],
				con = file.path(ZidDir, paste(VignNames[i], extension,
				sep = ".")))
		}
		# Rdata
		ZI.sample <- Zidb$.Data
		message("Extracting Rdata file...")
		save(ZI.sample, file = file.path(ZidDir, paste(sub(".zidb", "",
			basename(zidbfile)), "_dat1.RData", sep = "")))    
    
		# .zis data
		message("Extraction of ZIS data not supported yet...")
		## TODO...
	
		# Create zid file
		message("Compressing ZID file...")
		res <- zidCompress(zidir = ZidDir, type = "ZI3", check = FALSE,
			check.vignettes = FALSE, replace = replace, delete.source = TRUE)
	} else res <- TRUE
	
	# Do we have to delete the zidbfile?
	if (res && isTRUE(as.logical(delete.source))) unlink(zidbfile)
	
	message("-- Done! --")
		
	invisible(res)
}

# Convert .zidb files to .zid files
zidbToZidAll <- function (path = ".", zidbfiles, zisfiles =
file.path(dirname(zidbfiles), "Description.zis"), replace = FALSE,
delete.source = replace)
{
    ## First, switch to that directory
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- "."	# Indicate we are now in the right path
	
	## Get the list of zidbfiles to process
	if (missing(zidbfiles) || !length(zidbfiles))	# Compute them from path
		zidbfiles <- dir(path, pattern = extensionPattern("zidb"),
			full.names = TRUE) # All .zidb files

	## If there is no zidbfiles in this dir, exit now
	if (!length(zidbfiles)) {
		warning("There is no ZIDB files to process in ", getwd())
		return(invisible(FALSE))	
	}

	## Make sure there is no path associated
	#if (!all(zidbfiles == basename(zidbfiles))) {
	#	warning("You cannot provide paths for .zidb files, just file names")
	#	return(invisible(FALSE))
	#}
	
	## Check zisfiles and make sure the vector has same length as zidbfles
	## possibly recycling the file(s)
	zisfiles <- as.character(zisfiles)
	if (!length(zisfiles)) {
		warning("You must provide at least one ZIS file with samples characteristics")
		return(invisible(FALSE))
	}
	if (!checkFileExists(zisfiles, "zis", force.file = TRUE))
		return(invisible(FALSE))
	zisfiles <- rep(zisfiles, length.out = length(zidbfiles))
	
	## Create the .zid files from the .zidb files
	message("Conversion of ZIDB to ZID files...")
	flush.console()
	zidConvertOne <- function (item, zidbfiles, zisfiles, replace, delete.source)
		zidbToZid(zidbfiles[item], zisfile = zisfiles[item], replace = replace,
			delete.source = delete.source)
	items <- 1:length(zidbfiles)
	ok <- batch(items, zidConvertOne,zidbfiles = zidbfiles, zisfiles = zisfiles,
		replace = replace, delete.source = delete.source, verbose = FALSE)
	if (!ok) {
		warning(sum(attr(ok, "ok")), "/", length(zidbfiles),
			" items were correctly converted (see .last.batch)")
		invisible(FALSE)
	} else {
		message("-- Done! --")
		invisible(TRUE)
	}
}

## Link the database to R objects
zidbLink <- function (zidbfile)
	db2env(dbInit(zidbfile))

## Read only Rdata file from a .zidb database
zidbDatRead <- function (zidbfile)
	zidbLink(zidbfile)$.Data

## Read only the sample data
zidbSampleRead <- function (zidbfile)
	zidbLink(zidbfile)$.SampleData

## Functions to plot a collage
zidbPlotNew <- function (main = "ZooImage collage", ...)
{
	par(mfrow = c(1, 1), mar = c(0.1, 0.1, 2.1, 0.1))
	plot(0:1, 0:1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
		xaxs = "i", yaxs = "i", xlim = 0:1, ylim = 0:1, bty = "o",
		main = main, ...)
}

## Function to get a vignette from the database, rescale it and draw it in its
## corresponding vignette item
zidbDrawVignette <- function (rawimg, type = "jpg", item, nx = 5, ny = 5,
vmar = 0.01)
{
	## Centers for each vignette, on a graph area of [0, 1] on x and y
	nv <- nx * ny
	## Coordinates for centers of each vignette area
	xc <- (1:nx) / nx - 1/(nx * 2)
	yc <- (ny:1) / ny - 1/(ny * 2) # Because we want to start at the top and it
	## is the higher coord x and y coordinates for each vignette (fill from left
	## to right and top to bottom)
	vcoord <- expand.grid(list(x = xc, y = yc))
	## Half width and half height of a vignette area
	vhw <- ((xc[2] - xc[1]) - vmar) / 2
	vhh <- ((yc[1] - yc[2]) - vmar) / 2

	## Coordinates of top-left and bottom-right for vignettes areas
	vtl <- vcoord
	vtl$x <- vtl$x - vhw
	vtl$y <- vtl$y + vhh
	vbr <- vcoord
	vbr$x <- vbr$x + vhw
	vbr$y <- vbr$y - vhh

	## rawimg is a raw object containing JPEG or PNG data
	## item is the number of vignette area where to draw the vignette
	item <- as.integer(item[1])
	if (item < 1 || item > length(vtl$x)) stop("Wrong vignette item number")

	## Conversion from a raw object to a displayable image is done using
	## readPNG() or readJPEG() from the png/jpeg packages... For fast
	## processing, use native format, but 16bit not accepted for PNG and there
	## is a problem in case of transparency channel (if any) in PNG images on
	## windows devices
	if (type == "png") {
		vigimg <- readPNG(rawimg, native = TRUE)
	} else vigimg <- readJPEG(rawimg, native = TRUE)
	vigdim <- dim(vigimg) # Dimensions of the image in pixels
	## Determine top-left and bottom-right points of vignette bounding rectangle
	## for optimum display...
	## top-left point is always from the grid
	xleft <- vtl$x[item]
	ytop <- vtl$y[item]

	## Size of internal collage area (which is [0,1] both in x and y) in pixels
	totpx <- dev.size(units = "px")
	plt <- par("plt")
	totpx[1] <- totpx[1] * (plt[2] - plt[1]) # Width of collage area in pixels
	totpx[2] <- totpx[2] * (plt[4] - plt[3]) # Height of collage are in pixels

	## Size of vignette areas in pixels
	vwpx <- vhw * 2 * totpx[1]
	vhpx <- vhh * 2 * totpx[2]

	## If the vignette is smaller than the area, it fits without rescaling!
	if (vigdim[2] <= vwpx && vigdim[1] <= vhpx) {
		xright <- xleft + 2 * vhw / vwpx * vigdim[2]
		ybottom <- ytop - 2 * vhh / vhpx * vigdim[1]
	} else { # We need to rescale down the vignette to fit it in the area
		## Which dimension will fit the whole area?
		vigratio <- vigdim[2] / vigdim[1]
		arearatio <- vwpx / vhpx
		if (vigratio < arearatio) { # Fit height
			ybottom <- ytop - (2 * vhh)
			xright <- xleft + (2 * vhh * vigratio / arearatio)
		} else { # Fit width
			xright <- xleft + (2 * vhw)
			ybottom <- ytop - (2 * vhw / vigratio * arearatio)
		}
	}

	## Interpolation only works outside of windows!
	interpolate <- (names(dev.cur()) != "windows")

	## Note that if there is a transparency layer, a special treatment
	## is required for windows devices, see ?readPNG

	## Now, display that vignette in the collage
	rasterImage(vigimg, xleft, ybottom, xright, ytop, interpolate = interpolate)
}
