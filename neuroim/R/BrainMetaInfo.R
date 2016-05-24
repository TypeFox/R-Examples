

#' @include AllClass.R
roxygen()
#' @include AllGeneric.R
roxygen()
#' @include BrainFileDescriptor.R
roxygen()
#' @include Axis.R
roxygen()
#' @include NIFTI_IO.R
roxygen()

#' Generic function to create data reader
#' @param x an object specifying the infromation required to produce the reader
#' @param offset the byte offset (number of bytes to skip before reading)
#' @export dataReader
#' @rdname dataReader-methods
setGeneric(name="dataReader", def=function(x, offset) standardGeneric("dataReader"))

#' dim of \code{FileMetaInfo}
#' @param x the object
#' @export
setMethod(f="dim", signature=signature("FileMetaInfo"), 
		def=function(x) {
			x@Dim
		})


# @rdname loadData-methods
# setMethod(f="loadData", signature=signature(""))


#' @rdname dataReader-methods
setMethod(f="dataReader", signature=signature("NIfTIMetaInfo"), 
		def=function(x, offset=0) {
			if (x@fileDescriptor@dataEncoding == "gzip") {
				BinaryReader(gzfile(x@dataFile, "rb"), x@dataOffset+offset, .getRStorage(x@dataType), x@bytesPerElement, x@endian)	
			} else {
				BinaryReader(x@dataFile, x@dataOffset+offset, .getRStorage(x@dataType), x@bytesPerElement, x@endian)
			}
		})

#' @rdname dataReader-methods
setMethod(f="dataReader", signature=signature("AFNIMetaInfo"), 
		def=function(x, offset=0) {
			if (x@fileDescriptor@dataEncoding == "gzip") {
				BinaryReader(gzfile(x@dataFile, "rb"), x@dataOffset+offset, .getRStorage(x@dataType), x@bytesPerElement, x@endian)	
			} else {
				BinaryReader(x@dataFile, x@dataOffset+offset, .getRStorage(x@dataType), x@bytesPerElement, x@endian)
			}
		})		

#' @rdname dataReader-methods
setMethod(f="dataReader", signature=signature("NIMLSurfaceDataMetaInfo"), 
          def=function(x) {
            reader <- function(i) {
              if (length(i) == 1 && i == 0) {
                x@nodeIndices
              } else {
                x@data[,i,drop=FALSE]
              }
            }

            new("ColumnReader", nrow=as.integer(nrow(x@data)), ncol=as.integer(ncol(x@data)), reader=reader)
          })


#' @rdname readColumns-methods
setMethod(f="readColumns", signature=signature(x="ColumnReader", columnIndices="numeric"),
          def=function(x,columnIndices) {
            x@reader(columnIndices)
          })
          
            
#' @rdname trans-methods
setMethod(f="trans", signature=signature("BrainMetaInfo"), 
		def=function(x) {
			D <- min(length(x@Dim), 3)
			trans <- diag(c(x@spacing,1))
			trans[1:D,D+1] <- x@origin	
			trans
		})

#' @rdname trans-methods
setMethod(f="trans", signature=signature("NIfTIMetaInfo"), 
		def=function(x) {
			x@header$qform
		})

niftiDim <- function(nifti_header) {
	dimarray <- nifti_header$dimensions
	lastidx <- min(which(dimarray == 1)) - 1
	dimarray[2:lastidx]
}



#' This class contains meta information for an image
#'
#' @param Dim image dimensions
#' @param spacing voxel dimensions
#' @param origin coordinate origin
#' @param dataType the type of the data (e.g. "FLOAT")
#' @param label name(s) of images 
#' @param spatialAxes image axes for spatial dimensions (x,y,z)
#' @param additionalAxes axes for dimensions > 3 (e.g. time, color band, direction)
#' @return an instance of class \code{\linkS4class{BrainMetaInfo}}
#' @export BrainMetaInfo
#' @rdname BrainMetaInfo-class
BrainMetaInfo <- function(Dim, spacing, origin=rep(0, length(spacing)), dataType="FLOAT", label="", spatialAxes=OrientationList3D$AXIAL_LPI, additionalAxes=NullAxis) {
	new("BrainMetaInfo",
			Dim=Dim,
			spacing=spacing,
			origin=origin,
			dataType=dataType,
			label=label,
			spatialAxes=spatialAxes,
			additionalAxes=additionalAxes)
}						

#' Constructor for \code{\linkS4class{SurfaceGeometryMetaInfo}} class
#' @param descriptor the file descriptor
#' @param header a \code{list} containing header information
SurfaceGeometryMetaInfo <- function(descriptor, header) {
  stopifnot(is.numeric(header$vertices))
  stopifnot(is.numeric(header$faces))
 
  
  new("SurfaceGeometryMetaInfo",
     headerFile=header$headerFile,
     dataFile=header$dataFile,
     fileDescriptor=descriptor,
     vertices=as.integer(header$vertices),
     faces=as.integer(header$faces),
     label=as.character(header$label),
     embedDimension=as.integer(header$embedDimension))
}

#' Constructor for \code{\linkS4class{SurfaceDataMetaInfo}} class
#' @param descriptor the file descriptor
#' @param header a \code{list} containing header information
SurfaceDataMetaInfo <- function(descriptor, header) {
  stopifnot(is.numeric(header$nodes))
 
  new("SurfaceDataMetaInfo",
      headerFile=header$headerFile,
      dataFile=header$dataFile,
      fileDescriptor=descriptor,
      nodeCount=as.integer(header$nodes),
      nels=as.integer(header$nels),
      label=as.character(header$label))
}




#' Constructor for \code{\linkS4class{NIMLSurfaceDataMetaInfo}} class
#' @param descriptor the file descriptor
#' @param header a \code{list} containing header information
#' 
NIMLSurfaceDataMetaInfo <- function(descriptor, header) {
  stopifnot(is.numeric(header$nodes))
  
  new("NIMLSurfaceDataMetaInfo",
      headerFile=header$headerFile,
      dataFile=header$dateFile,
      fileDescriptor=descriptor,
      nodeCount=as.integer(header$nodeCount),
      nels=as.integer(header$nels),
      label=as.character(header$label),
      data=header$data,
      nodeIndices=header$nodes)
}

#' Constructor for \code{\linkS4class{NIfTIMetaInfo}} class
#' @param descriptor an instance of class \code{\linkS4class{NIfTIFileDescriptor}}
#' @param nifti_header a \code{list} returned by \code{readNIftiHeader}
#' @return an instance of class \code{\linkS4class{NIfTIMetaInfo}}
#' @export NIfTIMetaInfo
#' @rdname NIfTIMetaInfo-class
NIfTIMetaInfo <- function(descriptor, nifti_header) {
	stopifnot(!is.null(nifti_header$fileType) || (nifti_header$fileType == "NIfTI"))
	

	new("NIfTIMetaInfo",
			headerFile=headerFile(descriptor, nifti_header$fileName),
			dataFile=dataFile(descriptor, nifti_header$fileName),
			fileDescriptor=descriptor,
			endian=nifti_header$endian,
			dataOffset=nifti_header$voxOffset,
			dataType=nifti_header$dataStorage,
			bytesPerElement=as.integer(.getDataSize(nifti_header$dataStorage)),
			Dim=niftiDim(nifti_header),
			spatialAxes=.nearestAnatomy(nifti_header$qform),
			additionalAxes=NullAxis,
			spacing=nifti_header$pixdim[2:4],
			origin=nifti_header$qoffset,
			label=stripExtension(descriptor, basename(nifti_header$fileName)),
			intercept=nifti_header$sclIntercept,
			slope=nifti_header$sclSlope,
			header=nifti_header)
}

#' show an \code{SurfaceGeometryMetaInfo}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("SurfaceGeometryMetaInfo"), 
          def=function(object) {
            cat("an instance of class",  class(object), "\n\n")
            cat("number of vertices:", "\t", object@vertices, "\n")
            cat("number of faces:", "\t", object@faces, "\n")
            cat("label:", "\t", object@label, "\n")
            cat("embed dimension:", "\t", object@embedDimension, "\n")
          })

#' show an \code{SurfaceDataMetaInfo}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("SurfaceDataMetaInfo"), 
          def=function(object) {
            cat("an instance of class",  class(object), "\n\n")
            cat("nodeCount:", "\t", object@nodeCount, "\n")
            cat("nels:", "\t", object@nels, "\n")
            cat("label:", "\t", object@label, "\n")
          })

#' show a \code{FileMetaInfo}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("FileMetaInfo"), 
		def=function(object) {
			cat("an instance of class",  class(object), "\n\n")
			cat("headerFile:", "\t", object@headerFile, "\n")
			cat("dataFile:", "\t", object@dataFile, "\n")
			cat("endian:", "\t", object@endian, "\n")
			cat("dataOffset:", "\t", object@dataOffset, "\n")
			cat("dataType:", "\t", object@dataType, "\n")
			cat("dimensions:", "\t", object@Dim, "\n")
			cat("voxel size:", "\t", object@spacing, "\n")
			cat("origin:", "\t", object@origin, "\n")
			cat("label(s):", "\t", object@label, "\n")
			cat("intercept:", "\t", object@intercept, "\n")
			cat("slope:", "\t\t", object@slope, "\n\n")
			
			cat("additional format-specific info may be contained in @header slot", "\n")			
		})


#' AFNIMetaInfo 
#' 
#' Constructor for \code{\linkS4class{AFNIMetaInfo}} class
#' @param descriptor an instance of class \code{\linkS4class{AFNIFileDescriptor}}
#' @param afni_header a \code{list} returned by \code{readAFNIHeader}
#' @return an instance of class \code{\linkS4class{AFNIMetaInfo}}
#' @export AFNIMetaInfo
#' @rdname AFNIMetaInfo-class
AFNIMetaInfo <- function(descriptor, afni_header) {
		.Dim <- afni_header$DATASET_DIMENSIONS$content[afni_header$DATASET_DIMENSIONS$content > 0]
		if (afni_header$DATASET_RANK$content[2] > 1) {
			.Dim <- c(.Dim, afni_header$DATASET_RANK$content[2])			
		}
		
		
		labs <- if (is.null(afni_header$BRICK_LABS$content)) {
			labs <- paste("#", seq(0, afni_header$DATASET_RANK$content[2]), sep="")
		} else {
			afni_header$BRICK_LABS$content
		}
    
    
    ## AFNI contains a transform from IJK to dicom (RAI) space.
    ## We want the transform to go from IJK to nifti (LPI) space
		Tdicom <- matrix(afni_header$IJK_TO_DICOM$content, 3,4, byrow=TRUE)
    Tdicom <- rbind(Tdicom, c(0,0,0,1))
		TLPI <- diag(c(-1,-1,1,1)) %*% Tdicom
		
		new("AFNIMetaInfo",
			headerFile=headerFile(descriptor, afni_header$fileName),
			dataFile=dataFile(descriptor, afni_header$fileName),
			fileDescriptor=descriptor,
			endian=ifelse(afni_header[["BYTEORDER_STRING"]]$content == "MSB_FIRST", "big", "little"),
			dataOffset=0,
			dataType=switch(afni_header$BRICK_TYPES$content[1], "0"="BYTE", "1"="SHORT", "3"="FLOAT"),
			bytesPerElement=as.integer(switch(afni_header$BRICK_TYPES$content[1], "0"=1, "1"=2, "3"=4)),
			Dim=.Dim,
			spatialAxes=OrientationList3D$AXIAL_LPI,   # incorrect
			additionalAxes=NullAxis,            # incorrect
			spacing=abs(afni_header$DELTA$content),
			origin=afni_header$ORIGIN$content,
			label=labs,
			intercept=0,
			slope=ifelse(afni_header$BRICK_FLOAT_FACS$content == 0, 1, afni_header$BRICK_FLOAT_FACS$content),
			header=afni_header)
}
			
			
#' read header information of an image file
#'
#'
#' @param fileName the name of the file to read
#' @return an instance of class \code{\linkS4class{FileMetaInfo}} 
#' @export readHeader
readHeader <- function(fileName) {
	desc <- findDescriptor(fileName) 
	if (is.null(desc)) {
		stop(paste("could not find reader for file: ", fileName))
	}
	
	readMetaInfo(desc, fileName)			
}

setAs(from="BrainMetaInfo", to="NIfTIMetaInfo", def=function(from) {
			if (inherits(from, "NIfTIMetaInfo")) {
				from
			} else {
				hdr <- as.nifti.header(from)
				desc <- findDescriptor(hdr$fileName)
				NIfTIMetaInfo(desc, hdr)
			}
						
		})



