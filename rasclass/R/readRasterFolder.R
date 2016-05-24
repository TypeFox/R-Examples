##################################################################################
# Set Method: readRasterFolder
##################################################################################
readRasterFolder <- function(path, samplename = 'sample', filenames = NULL, object = new('rasclass'), asInteger = FALSE){}

setMethod('readRasterFolder', signature(path = 'character'),

function(path, samplename = 'sample', filenames = NULL, object = new('rasclass'),  asInteger = FALSE){

	# Clean and store samplename
	if(substr(samplename, nchar(samplename)-3, nchar(samplename)) == '.asc'){
		samplename <- substr(samplename, 1, nchar(samplename)-4)
	}
	object@samplename <- samplename

	# Clean and store path
	if(substr(path,nchar(path), nchar(path)) != '/'){
		path <- paste(path,'/', sep='')
	}
	object@path <- path

	# Remember current working directory and set it back on exit, then
	# change to data folder.
	userwd <- getwd()
	on.exit(setwd(userwd))
	setwd(path)

	# Read Variable Names
	if(length(filenames) == 0){
		filelist <- Sys.glob('*.asc')
		filelist <- filelist[filelist != paste(object@samplename, '.asc' , sep='')]
	}
	else{
		filelist <- NA
		for(i in 1:length(filenames)){
			file <- filenames[i]
			if(substr(file, nchar(file)-3, nchar(file)) == '.asc'){
				filelist[i] <- file
			} else {
				filelist[i] <- paste(file, '.asc', sep='')
			}
		}
	}

	# Store file names as names for data columns
	namelist <- substr(filelist, 1, nchar(filelist)-4)

	# Load Sample
	cat('\nReading Raster grids from "', getwd(), '"\n', sep ='')
	cat(paste(paste(object@samplename),'.asc',sep=''), '\n')

	thissample  <- readRaster(paste(object@samplename, '.asc', sep=''), asInteger = asInteger)
	object@data <- data.frame(thissample@grid)

	# Set gridSkeleton headers
	object@gridSkeleton@ncols 		<- thissample@ncols
	object@gridSkeleton@nrows 		<- thissample@nrows
	object@gridSkeleton@xllcorner 	<- thissample@xllcorner
	object@gridSkeleton@yllcorner 	<- thissample@yllcorner
	object@gridSkeleton@cellsize 	<- thissample@cellsize
	object@gridSkeleton@NAvalue 	<- thissample@NAvalue
	rm(thissample)

	# Load Dependent Varialbe Data using Filenames
	for(i in 1:length(filelist)){

		# Communicate file to read
		cat(filelist[i], '\n')

		# Load data from file
		tempraster <- readRaster(filelist[i], asInteger = asInteger)

		# Check if header is equal to sample
		if( object@gridSkeleton@nrows     != tempraster@nrows |
			object@gridSkeleton@ncols     != tempraster@ncols |
			object@gridSkeleton@xllcorner != tempraster@xllcorner |
			object@gridSkeleton@yllcorner != tempraster@yllcorner |
			object@gridSkeleton@cellsize  != tempraster@cellsize
			){ stop('The input raster grids do not all have the same header') }

		# Update gridSkeleton and drop uncomplete rows
		if(i == 1){
			# Add data to dataframe
			object@data <- data.frame(object@data, tempraster@grid)

			# Create grid skeleton
			object@gridSkeleton@grid <- as.integer(!is.na(tempraster@grid))
		}
		else{
			# Add data to dataframe, omitting unused rows with old gridSkeleton
			object@data <- data.frame(object@data, tempraster@grid[as.logical(object@gridSkeleton@grid)])

			# Update gridSkeleton
			object@gridSkeleton@grid <- as.integer(object@gridSkeleton@grid * !is.na(tempraster@grid))
		}

		# Remove rows with NAs in the new column
		object@data <- object@data[!is.na(object@data[, i+1]), ]
	}

	# Set Names in Data Frame
	names(object@data) <- append(object@samplename, namelist)

	# Build Formula
	object <- buildFormula(object)

	# Return object
	object
}
)
