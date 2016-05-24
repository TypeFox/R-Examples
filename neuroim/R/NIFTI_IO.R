#' @include common.R
{}
#' @include BinaryIO.R
{}


.checkDimensions <- function(dimvec) {
	if (any(dimvec < 0)) {
		stop(paste("nifti(checkDimensons): illegal dimension vector in header: ", dimvec))
	}
}

write.nifti.vector <- function(vec, fileName, dataType=NULL) {      
	stopifnot(length(dim(vec)) == 4)
	hdr <- as.nifti.header(vec, vec@source@metaInfo)
	
	if (!is.null(dataType) && dataType != vec@source@metaInfo@dataType) {
		hdr$datatype <- .getDataCode(dataType)
		hdr$dataStorage <- .getDataStorage(hdr$datatype)		
		hdr$bitpix <- .getDataSize(dataType) * 8		
	} else {
		dataType <- vec@source@metaInfo@dataType
		
		### code duplication
		hdr$datatype <- .getDataCode(dataType)
		hdr$dataStorage <- .getDataStorage(hdr$datatype)		
		hdr$bitpix <- .getDataSize(dataType) * 8		
	}
	
	conn <- if (substr(fileName, nchar(fileName)-2, nchar(fileName)) == ".gz") {
				gzfile(fileName, open="wb")
			} else {
				file(fileName, open="wb")
			}
	
	writeNIfTIHeader(hdr, conn, close=FALSE)
	writer <- BinaryWriter(conn, hdr$voxOffset, dataType, hdr$bitpix/8, .Platform$endian)
	
	NVOLS <- dim(vec)[4]
	for (i in 1:NVOLS) {
		writeElements(writer, as.numeric(takeVolume(vec, i)))
	}
	close(writer)
}  


write.nifti.volume <- function(vol, fileName, dataType=NULL) {      
	stopifnot(length(dim(vol)) == 3)
	hdr <- as.nifti.header(vol, vol@source@metaInfo)
	
	if (!is.null(dataType) && dataType != vol@source@metaInfo@dataType) {
		hdr$datatype <- .getDataCode(dataType)
		hdr$dataStorage <- .getDataStorage(hdr$datatype)		
		hdr$bitpix <- .getDataSize(dataType) * 8		
	} else {
		dataType <- vol@source@metaInfo@dataType
		### code duplication
		hdr$datatype <- .getDataCode(dataType)
		hdr$dataStorage <- .getDataStorage(hdr$datatype)		
		hdr$bitpix <- .getDataSize(dataType) * 8		
	}
	
	conn <- if (substr(fileName, nchar(fileName)-2, nchar(fileName)) == ".gz") {
		gzfile(fileName, open="wb")
	} else {
		file(fileName, open="wb")
	}

	writeNIfTIHeader(hdr, conn, close=FALSE)
	writer <- BinaryWriter(conn, hdr$voxOffset, dataType, hdr$bitpix/8, .Platform$endian)
	writeElements(writer, as.numeric(vol))
	close(writer)
}  



as.nifti.header <- function(vol, metaInfo, fname=NULL, oneFile=TRUE) {
	if (inherits(metaInfo, "NIfTIMetaInfo")) {
		ret <- metaInfo@header
		if (!is.null(fname)) {
			ret$fileName <- fname
		}
		ret$onefile <- oneFile
		if (oneFile) {
			ret$magic <- "n+1"
		} else {
			ret$magic <- "ni1"
		}
		
		ret
	} else if (inherits(metaInfo, "FileMetaInfo")){
		### serious code duplication
		hd <- createNIfTIHeader(oneFile=oneFile, fileName=fname)
		if (is.null(fname)) {
			hd$fileName <- metaInfo@headerFile
		} else {
			hd$fileName <- fname
		}
	
		hd$endian <- metaInfo@endian
		hd$voxOffset <- metaInfo@dataOffset
		hd$datatype <- .getDataCode(metaInfo@dataType)
		hd$dataStorage <- .getDataStorage(hd$datatype)
		hd$bitpix <- metaInfo@bytesPerElement * 8
		hd$dimensions <- c(length(metaInfo@Dim), metaInfo@Dim)
		N <- 8 - length(hd$dimensions)
		hd$dimensions <- c(hd$dimensions,  rep(1, N))
		hd$numDimensions <- length(metaInfo@Dim)
		
		### only encodes pixdim for three dimensions
		hd$pixdim <- c(0, metaInfo@spacing, rep(0,4))
		hd$qoffset <- metaInfo@origin
		hd$sclIntercept <- metaInfo@intercept
		hd$sclSlope <- metaInfo@sclSlope
	
		#tmat <- diag(c(metaInfo@spacing, 1))
		#tmat[1:3, 4] <- metaInfo@origin
		
		
		
		tmat <- trans(vol)
		
		hd$qform <- tmat
		hd$sform <- tmat
		
		quat1 <- .matrixToQuatern(tmat)
		hd$quaternion <- quat1$quaternion
		hd$qfac <- quat1$qfac
		hd$pixdim[1] <- hd$qfac
		hd
		### serious code duplication
	} else if (inherits(metaInfo, "BrainMetaInfo")) {
		hd <- createNIfTIHeader(oneFile=oneFile, fileName=fname)
		
		if (is.null(fname)) {
			hd$fileName <- "nothing.nii"
		} else {
			hd$fileName <- fname
		}
		
		
		### serious code duplication
		
		hd$endian <- .Platform$endian
		hd$voxOffset <- 352
		
		hd$datatype <- .getDataCode(metaInfo@dataType)
		hd$dataStorage <- .getDataStorage(hd$datatype)
		
		hd$bitpix <- .getDataSize(metaInfo@dataType) * 8
				
				
		hd$dimensions <- c(length(metaInfo@Dim), metaInfo@Dim)
		N <- 8 - length(hd$dimensions)
		hd$dimensions <- c(hd$dimensions,  rep(1, N))
		hd$numDimensions <- length(metaInfo@Dim)
		
		### only encodes pixdim for three dimensions
		hd$pixdim <- c(0, metaInfo@spacing, rep(0,4))
		
		
		hd$qoffset <- metaInfo@origin
		hd$sclIntercept <- 0
		hd$sclSlope <- 0
		
		
		tmat <- trans(vol)
		
		
		hd$qform <- tmat
		hd$sform <- tmat
		
		quat1 <- .matrixToQuatern(tmat)
		hd$quaternion <- quat1$quaternion
		hd$qfac <- quat1$qfac
		hd$pixdim[1] <- hd$qfac
		hd
		
		### serious code duplication
		
		
	}
		
}

createNIfTIHeader <- function(oneFile=TRUE, fileName=NULL) {
	header <- list()
	header$fileType <- "NIfTI"
	header$encoding <- "binary"
	header$version <- "1"


	header$fileName <- fileName
	header$endian <- .Platform$endian
	
	header$diminfo <- 0
	header$dimensions <- NULL
	header$numDimensions <- NULL
	
	header$intent1 <-  0
	header$intent2 <-  0
	header$intent3 <-  0
	
	
	header$intentCode <-  0
	header$datatype <- NULL
	header$dataStorage <- NULL
	header$bitpix <- NULL
	header$sliceStart <- 0
	header$pixdim <-  NULL
	
	header$qfac <- -1
	
	
	header$voxOffset <- 0
	header$sclSlope <- 0
	header$sclIntercept <- 0 
	header$sliceEnd <- 0
	header$sliceCode <- 0
	
	
	header$xyztUnits <- 2
	header$calMax <- 0
	header$calMin <- 0
	
	header$sliceDuration <- 0
	header$toffset <- 0
	
	header$glmax <- 0
	header$glmin <- 0
	
	header$description <- character(80)
	header$auxfile <- character(24)
	
	header$qformCode <- 1
	header$sformCode <- 1
	header$quaternion <- NULL
	
	header$qoffset <- NULL
	header$qform <- NULL
	
	header$sform <- NULL
	header$intentName <- character(16)
	header$magic <- "n+1"
	
	header$onefile <- oneFile
	if (oneFile) {
		header$magic <- "n+1"
	} else {
		header$magic <- "ni1"
	}
	
	header$version <- 1
	header
	
}

readNIfTIHeader <- function(fname) {
	

	header <- list()
	header$fileType <- "NIfTI"
	header$encoding <- "binary"
	header$version <- "1"
	
	conn <- NULL
	
	if (.isExtension(fname, ".nii") || .isExtension(fname, ".hdr")) {
		conn <- file(fname, open="rb")
	} else if (.isExtension(fname, ".nii.gz")) {
		conn <- gzfile(fname, open="rb")
		header$encoding <- "gzip"
	} else {
		stop(paste("illegal NIFTI header name", fname))
	}
	
	
	endian <- .getEndian(conn)
	
	header$fileName <- fname
	header$endian <- endian
	
	readBin(conn, what=integer(), n=10+18+4+2+1, size=1)
	
	header$diminfo <- readBin(conn, what=integer(), n=1, size=1)
	header$dimensions <- readBin(conn, integer(), n=8, size=2, endian=endian)
	
	header$dimensions[header$dimensions == 0] <- 1
	
	.checkDimensions(header$dimensions)
	
	
	header$numDimensions <- header$dimensions[1]
	
	header$intent1 <-  readBin(conn, double(), n=1, size=4, endian=endian)
	header$intent2 <-  readBin(conn, double(), n=1, size=4, endian=endian)
	header$intent3 <-  readBin(conn, double(), n=1, size=4, endian=endian)
	
	
	header$intentCode <-  readBin(conn, integer(), n=1, size=2, endian=endian)
	header$datatype <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$dataStorage <- .getDataStorage(header$datatype)
	header$bitpix <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$sliceStart <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$pixdim <-  readBin(conn, double(), n=8, size=4, endian=endian)

	header$qfac = header$pixdim[1]
	
	if (header$qfac == 0) {
		header$qfac = 1
	}
	
	header$voxOffset <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$sclSlope <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$sclIntercept <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$sliceEnd <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$sliceCode <-  readBin(conn, integer(), n=1, size=1, endian=endian)
	
	
	header$xyztUnits <- readBin(conn, integer(), n=1, size=1, endian=endian)
	header$calMax <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$calMin <- readBin(conn, double(), n=1, size=4, endian=endian)
	
	header$sliceDuration <- readBin(conn, double(), n=1, size=4, endian=endian)
	header$toffset <- readBin(conn, double(), n=1, size=4, endian=endian)
	
	header$glmax <- readBin(conn, integer(), n=1, size=4, endian=endian) # unused glmax, glmin
	header$glmin <- readBin(conn, integer(), n=1, size=4, endian=endian) # unused glmax, glmin
	
	header$description <- readBin(conn, integer(), n=80, size=1, endian=endian)
	header$auxfile <- readBin(conn, integer(), n=24, size=1, endian=endian)
	
	header$qformCode <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$sformCode <- readBin(conn, integer(), n=1, size=2, endian=endian)
	header$quaternion <- readBin(conn, double(), n=3, size=4, endian=endian)
	
	header$qoffset <- readBin(conn, double(), n=3, size=4, endian=endian)
	header$qform <- .quaternToMatrix(header$quaternion, header$qoffset, header$pixdim[2:4], header$qfac)
	
	sform  <- readBin(conn, double(), n=12, size=4, endian=endian)
	header$sform <- rbind(matrix(sform,3,4, byrow=T), c(0,0,0,1))
	header$intentName <- readBin(conn, character(), n=16, size=1, endian=endian)
	header$magic <- readChar(conn, nchars=4)
	
	header$onefile <- F
	if (substr(header$magic,2,2) == "+") {
		header$onefile <- T
	}
	
	header$version <- substr(header$magic,3,3)
	
	close(conn)
	
	header
	
}

writeNIfTIHeader <- function(niftiInfo, conn, close=TRUE) {
	endian <- niftiInfo$endian
	
		
	writeBin(as.integer(348), conn, 4, endian) 
	writeBin(integer(34),conn,1,endian)
	writeChar("r", conn,1,eos=NULL)
	writeBin(as.integer(niftiInfo$diminfo), conn, size=1, endian) #diminfo, not supported currently -- write zero  
	#writeBin(as.integer(niftiInfo$numDimensions), conn, 2, endian)         #num dimensions 
	
	stopifnot(length(niftiInfo$dimensions) == 8)
	stopifnot(length(niftiInfo$pixdim) == 8)
	
	writeBin(as.integer(niftiInfo$dimensions), conn, 2, endian)   #dimension vector 
	writeBin(as.double(niftiInfo$intent1), conn, 4, endian)       #intent1
	writeBin(as.double(niftiInfo$intent2), conn, 4, endian)       #intent2
	writeBin(as.double(niftiInfo$intent3), conn, 4, endian)       #intent3
	writeBin(as.integer(niftiInfo$intentCode), conn, 2, endian)   #intent code
	writeBin(as.integer(niftiInfo$datatype),conn, 2, endian)      #datatype 
	writeBin(as.integer(niftiInfo$bitpix),conn, 2, endian)        #bits per pixel 
	writeBin(as.integer(niftiInfo$sliceStart),conn, 2, endian)    #slice start
	writeBin(as.double(niftiInfo$pixdim), conn, 4, endian)        #pix dim
	writeBin(as.double(niftiInfo$voxOffset), conn, 4, endian)     #voxel offset
	writeBin(as.double(niftiInfo$sclSlope), conn, 4, endian)      #slope
	writeBin(as.double(niftiInfo$sclIntercept), conn, 4, endian)  #intercept
	writeBin(as.integer(niftiInfo$sliceEnd), conn, 2, endian)     #slice end
	writeBin(as.integer(niftiInfo$sliceCode), conn, 1, endian)    #slice code
	writeBin(as.integer(niftiInfo$xyztUnits), conn, 1, endian)    #xyzt units
	writeBin(as.double(niftiInfo$calMax), conn, 4, endian)        #cal max
	writeBin(as.double(niftiInfo$calMin), conn, 4, endian)        #cal min
	writeBin(as.double(niftiInfo$sliceDuration), conn, 4, endian) #slice duration
	writeBin(as.double(niftiInfo$toffset), conn, 4, endian)       #t offset
	writeBin(as.integer(niftiInfo$glmax), conn, 4, endian)        #glmax
	writeBin(as.integer(niftiInfo$glmin), conn, 4, endian)        #glmin
	writeBin(as.integer(niftiInfo$description), conn, 1, endian)  #description
	writeBin(as.integer(niftiInfo$auxfile), conn, 1, endian)      #aux_file
	writeBin(as.integer(niftiInfo$qformCode), conn, 2, endian)    #qform code
	writeBin(as.integer(niftiInfo$sformCode), conn, 2, endian)    #sform code
	
	writeBin(as.double(niftiInfo$quaternion), conn, 4, endian)    #quaternion
	writeBin(as.double(niftiInfo$qoffset), conn, 4, endian)       #qoffset
	writeBin(as.double(t(niftiInfo$sform[1:3,])), conn, 4, endian) #sform
	writeBin(as.integer(niftiInfo$intentName), conn, 1, endian)    #intentName
	writeChar(niftiInfo$magic, conn)                               #magic
	
	loc <- seek(conn)
	offset <- niftiInfo$voxOffset
	
	nbytes <- offset-loc
	
	## doesn't support extensions yet
	if (nbytes > 0) {
		writeBin(integer(nbytes), conn, size=1, endian)
	}
	
	if (close) {
		close(conn)
	}
	
	conn
}
