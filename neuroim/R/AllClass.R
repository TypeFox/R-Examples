#' @include AllGeneric.R
#' @import methods
roxygen()


setOldClass("mesh3d")
setOldClass(c("file", "connection"))
setOldClass(c("gzfile", "connection"))
setOldClass("raster")
setOldClass("rastergrob")
setOldClass("gList")

.package_env <- new.env()

#' Base
#' 
#' Generic S4 Base class
#' 
#' @rdname Base-class
setClass("Base", contains=c("VIRTUAL"))


#' NamedAxis
#' 
#' This class represents an axis with a name attribute
#' @rdname NamedAxis-class
#' @slot axis the name of the axis
#' @slot direction of axis (-1,+1)
#' @export 
setClass("NamedAxis", representation=
				representation(axis="character", direction="numeric"))

#' Base
#' 
#' Virtual base class representing an ordered set of named axes.
#' @rdname AxisSet-class
#' @slot ndim the number of axes (or dimensions)
#' @export
setClass("AxisSet", representation=representation(ndim="integer"))

#' AxisSet1D
#' 
#' A one-dimensional axis set
#' @rdname AxisSet1D-class
#' @slot i the first axis
#' @export
setClass("AxisSet1D", representation=representation(i="NamedAxis"), contains=c("AxisSet"))

#' AxisSet2D
#' 
#' A two-dimensional axis set
#'
#' @rdname AxisSet2D-class
#' @slot j the second axis
#' @export
setClass("AxisSet2D", representation=representation(j="NamedAxis"), 
		contains=c("AxisSet1D"))

#' AxisSet3D
#' 
#' A three-dimensional axis set
#' 
#' @rdname AxisSet3D-class
#' @slot k the third axis
#' @export
setClass("AxisSet3D", representation=representation(k="NamedAxis"),
		contains=c("AxisSet2D"))

#' AxisSet4D
#' 
#' A four-dimensional axis set
#' @name AxisSet4D-class
#' @slot l the fourth axis
#' @export
setClass("AxisSet4D", representation=representation(l="NamedAxis"),
		contains=c("AxisSet3D"))

#' AxisSet5D
#' 
#' A five-dimensional axis set
#'
#' @name AxisSet5D-class
#' @slot m the fifth axis
#' @export
setClass("AxisSet5D", representation=representation(m="NamedAxis"),
		contains=c("AxisSet4D"))

#' BrainFileDescriptor
#' 
#' This class represents a neuroimaging file format
#'
#' @rdname BrainFileDescriptor-class
#' @slot fileFormat the name of the file format (e.g. NIfTI)
#' @slot headerEncoding the file encoding of the header file (e.g. 'raw' for binary, 'gzip' for gz compressed')
#' @slot headerExtension the file extension for the header file (e.g. 'nii' for NIfTI single files)
#' @slot dataEncoding the file encoding for the data file
#' @slot dataExtension the file extension for the data file (e.g. 'nii' for NIfTI single files)
#' @exportClass BrainFileDescriptor
setClass("BrainFileDescriptor",
		representation=
				representation(fileFormat="character",
						headerEncoding="character",
						headerExtension="character",
						dataEncoding="character",
						dataExtension="character"),
    )

#' NIfTIFileDescriptor
#'
#' This class supports the NIfTI file format
#' 
#' @rdname NIfTIFileDescriptor-class
#' @export
setClass("NIfTIFileDescriptor", contains=c("BrainFileDescriptor"))

#' AFNIFileDescriptor
#' 
#' This class supports the AFNI file format
#' @rdname AFNIFileDescriptor-class
#' @export
setClass("AFNIFileDescriptor", contains=c("BrainFileDescriptor"))


#' NIMLSurfaceFileDescriptor
#' 
#' This class supports the NIML file format for surface-based data
#' @rdname NIMLSurfaceFileDescriptor-class
#' @export
setClass("NIMLSurfaceFileDescriptor", contains=c("BrainFileDescriptor"))

#' FresurferAsciiSurfaceFileDescriptor
#' 
#' This class supports the FreesurferAsciiSurfaceFileDescriptor file format for surface geometry
#' @rdname FreesurferAsciiSurfaceFileDescriptor-class
#' @export
setClass("FreesurferAsciiSurfaceFileDescriptor", contains=c("BrainFileDescriptor"))


#' BaseMetaInfo
#' 
#' This is a base class to represent meta information
#' @rdname BaseMetaInfo-class
#' @exportClass BaseMetaInfo
setClass("BaseMetaInfo")

#' NullMetaInfo
#' 
#' This is class is used to denote the absense of meta information
#' @rdname NullMetaInfo-class
setClass("NullMetaInfo", contains=c("BaseMetaInfo"))


#' show a \code{BaseMetaInfo}
#' @param object the object
#' @export
setMethod(f="show",
		signature=signature(object="BaseMetaInfo"),
		def=function(object) {
			cat("an instance of class",  class(object), "\n\n")
		})

#' show a \code{NullMetaInfo}
#' @param object the object
#' @export
setMethod(f="show",
		signature=signature(object="NullMetaInfo"),
			def=function(object) {
				cat("an instance of class",  class(object), "\n\n")
				cat("meta info is null \n")
			})
	

#' SurfaceGeometryMetaInfo
#' 	
#' This class contains meta information for brain surface geometry
#' 
#' @rdname SurfaceGeometryMetaInfo-class
#' @slot headerFile name of the file containing meta information
#' @slot dataFile name of the file containing data
#' @slot fileDescriptor descriptor of image file format
#' @slot vertices the number of surface vertices
#' @slot faces the number of faces
#' @slot embedDimension the dimensionality of the embedding
#' @slot label a label indicating the type of surface (e.g. white, pial, inflated, flat, spherical)
#' @export	 							 
setClass("SurfaceGeometryMetaInfo",
         representation=
           representation(
             headerFile="character",
             dataFile="character",
             fileDescriptor="BrainFileDescriptor",
             vertices="integer",
             faces="integer",
             label="character",
             embedDimension="integer"),
         contains=c("BaseMetaInfo"))

#' FreeSurferSurfaceGeometryMetaInfo
#' 	
#' This class contains meta information for brain surface geometry
#' 
#' @rdname FreeSurferSurfaceGeometryMetaInfo-class
#' @export	 							 
setClass("FreesurferSurfaceGeometryMetaInfo", contains=c("SurfaceGeometryMetaInfo"))

#' SurfaceDataMetaInfo
#' 	
#' This class contains meta information for surface-based data (the values that map to a surface geometry)
#' 
#' @rdname SurfaceDataMetaInfo-class
#' @slot headerFile name of the file containing meta information
#' @slot dataFile name of the file containing data
#' @slot fileDescriptor descriptor of image file format
#' @slot nodeCount the number of nodes for which surface data exists
#' @slot nels the number of data vectors (typically the number of columns in the surface data matrix; nels = 1 for a single surface data set)
#' @slot label a label indicating the type of surface (e.g. white, pial, inflated, flat, spherical)
#' @export	 							 
setClass("SurfaceDataMetaInfo",
         representation=
           representation(
             headerFile="character",
             dataFile="character",
             fileDescriptor="BrainFileDescriptor",
             nodeCount="integer",
             nels="integer",
             label="character"),
         contains=c("BaseMetaInfo"))

#' NIMLSurfaceDataMetaInfo
#' 	
#' This class contains meta information for surface-based data for the NIML data format
#' 
#' @rdname NIMLSurfaceDataMetaInfo-class
#' @slot data the numeric data matrix of surface values (rows = nodes, columns=surface vectors)
#' @slot nodeIndices the indices of the nodes for mapping to associated surface geometry.
#' @export	 							 
setClass("NIMLSurfaceDataMetaInfo",
         representation=
           representation(
             data="matrix",
             nodeIndices="integer"),
         contains=c("SurfaceDataMetaInfo"))

#' BrainMetaInfo
#' 	
#' This class contains meta information from an image
#' @rdname BrainMetaInfo-class
#' @slot dataType the data type code, e.g. FLOAT
#' @slot Dim image dimensions
#' @slot spatialAxes image axes for spatial dimensions (x,y,z)
#' @slot additionalAxes axes for dimensions > 3 (e.g. time, color band, direction)
#' @slot spacing voxel dimensions
#' @slot origin coordinate origin
#' @slot label name(s) of images 
#' @export	 							 
setClass("BrainMetaInfo",
			representation=
					representation(
							dataType="character",
							Dim="numeric",
							spatialAxes="AxisSet3D",
							additionalAxes="AxisSet",
							spacing="numeric",
							origin="numeric",
							label="character"),
			
			#prototype=prototype(),
			contains=c("BaseMetaInfo"))

#' FileMetaInfo
#' 
#' This class contains meta information from an image data file
#'
#' @rdname FileMetaInfo-class
#' @slot headerFile name of the file containing meta information
#' @slot dataFile name of the file containing data
#' @slot fileDescriptor descriptor of image file format
#' @slot endian byte order of data ('little' or 'big')
#' @slot dataOffset the number of bytes preceding the start of image data in data file
#' @slot bytesPerElement number of bytes per element
#' @slot intercept constant value added to image -- multiple values allowed (must equal numer of sub-images)
#' @slot slope image multiplier -- multiple values allowed (must equal numer of sub-images)  
#' @slot header a list of format specific attributes
#' @export	 							 
setClass("FileMetaInfo",
		representation=
				representation(headerFile="character",
					   dataFile="character",
					   fileDescriptor="BrainFileDescriptor",
					   endian="character",
					   dataOffset="numeric",					   
					   bytesPerElement="integer",
					   intercept="numeric",
					   slope="numeric",
					   header="list"),
					   
	    #prototype=prototype(),
		contains=c("BrainMetaInfo"))

#' NifTIMetaInfo
#' 
#' This class contains meta information for a NIfTI image file
#' @rdname FileMetaInfo-class
#' @slot nifti_header a \code{list} of attributes specific to the NIfTI file format 
#' @export	
setClass("NIfTIMetaInfo",
    representation=(nifti_header="list"),
		contains=c("FileMetaInfo"))

#' AFNIMetaInfo
#' 
#' This class contains meta information for a AFNI image file
#' @rdname FileMetaInfo-class  
#' @slot afni_header a list of attributes specific to the AFNI file format 
#' @slot afni_header a \code{list} of attributes specific to the AFNI file format 
#' @export
setClass("AFNIMetaInfo",
    representation=(afni_header="list"),
		contains=c("FileMetaInfo"))

#' BaseSource
#' 
#' This is a base class to represent a data source
#' @rdname BaseSource-class
#' @slot metaInfo meta information for the data source
#' @exportClass BaseSource
setClass("BaseSource", representation=representation(metaInfo="BaseMetaInfo"))


NullSource <- function() {
	ex <- exists(".NullSource", envir=.package_env)
	if (!ex) {
		ret <- new("BaseSource", metaInfo=new("NullMetaInfo"))	
		assign(".NullSource", ret, envir=.package_env)
		ret
	} else {
		get(".NullSource", envir=.package_env)
	}
	
}

#' BrainSource
#' 
#' Base class for representing a data source for images. The purpose of this class is to provide a layer in between 
#' low level IO and image loading functionality.
#' @rdname BrainSource-class
#' @slot metaInfo meta information for the data source
#' @exportClass BrainSource
setClass("BrainSource", representation=
				representation(metaInfo="BrainMetaInfo"),
				contains=c("BaseSource"))

#' BrainFileSource
#' 		
#' Base class for representing a data source for images. The purpose of this class is to provide a layer in between 
#' low level IO and image loading functionality.
#' @rdname BrainFileSource-class
#' @slot metaInfo meta information for the data source
#' @exportClass BrainFileSource
setClass("BrainFileSource", representation=
				representation(metaInfo="FileMetaInfo"),
				contains=c("BrainSource"))

		
#' BrainVolume
#' BrainVolumeSource
#' 		
#' A class is used to produce a \code{\linkS4class{BrainVolume}} instance
#' @rdname BrainVolumeSource-class
#' @slot index the index of the volume to be read -- must be of length 1.
#' @exportClass BrainVolumeSource
		setClass("BrainVolumeSource", representation=
						representation(index="integer"),
				contains=c("BrainSource"))

#' BrainVectorSource
#' 
#' A class that is used to produce a \code{\linkS4class{BrainVectorSource}} instance
#' @rdname BrainVectorSource-class
#' @slot indices the index vector of the volumes to be loaded
#' @export
		setClass("BrainVectorSource", representation=
						representation(indices="integer"),
				contains=c("BrainSource"))


#' BrainSurfaceVectorSource
#' 
#' A class that is used to produce a \code{\linkS4class{BrainSurfaceVectorSource}} instance
#' 
#' @rdname BrainSurfaceVectorSource-class
#' @slot indices the index vector of the volumes to be loaded
#' @export
setClass("BrainSurfaceVectorSource", representation=
           representation(indices="integer"),
         contains=c("BaseSource"))

		
		
#' BrainBucketSource
#' 
#' A class that is used to produce a \code{\linkS4class{BrainBucket}} instance
#' @rdname BrainBucketSource-class
#' @slot sourceList a list of sources for the bucket sub-volumes
#' @slot cache a cache used to store data in memory
#' @export
		setClass("BrainBucketSource",
				representation=representation(sourceList="list", cache="environment"),
				contains=c("BrainVectorSource"))
		
#' BrainSurfaceSource
#' 
#' A class that is used to produce a \code{\linkS4class{BrainSurface}} instance
#' @rdname BrainSurfaceSource-class
#' @slot metaInfo a \code{\linkS4class{SurfaceGeometryMetaInfo}} instance
#' @slot dataMetaInfo a \code{\linkS4class{SurfaceDataMetaInfo}} instance
#' @slot index the index offset into the surface data matrix
#' @export
setClass("BrainSurfaceSource",
         representation=representation(metaInfo="SurfaceGeometryMetaInfo", dataMetaInfo="SurfaceDataMetaInfo", index="integer"),
         contains=c("BaseSource"))


#' BinaryReader
#' 
#' This class supports reading of bulk binary data from a connection
#' @rdname BinaryReader-class
#' @slot input the binary input connection
#' @slot byteOffset the number of bytes to skip at the start of input
#' @slot dataType the dataType of the binary Elements
#' @slot bytesPerElement number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @slot endian endianness of binary input connection
#' @export
setClass("BinaryReader", representation=
				representation(input="connection",
							   byteOffset="numeric",
							   dataType="character",
							   bytesPerElement="integer",
							   endian="character"))


#' ColumnReader
#' 
#' This class supports reading of data froma matrix-like stroage format
#' @rdname ColumnReader-class
#' @slot nrow the number of rows
#' @slot ncol the number of columns
#' @slot reader a function that takes a set of column indices and returns a \code{matrix}
#' @export
setClass("ColumnReader", representation=
           representation(nrow="integer", ncol="integer", reader="function"))

         
#' BinaryWriter
#' 
#' This class supports writing of bulk binary data to a connection
#' @rdname BinaryWriter-class
#' @slot output the binary output connection
#' @slot byteOffset the number of bytes to skip at the start of input
#' @slot dataType the dataType of the binary Elements
#' @slot bytesPerElement number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @slot endian endianness of binary output connection
#' @export
setClass("BinaryWriter", representation=
							   representation(output="connection",	
							   byteOffset="numeric",
							   dataType="character",
							   bytesPerElement="integer",
							   endian="character"))
					   
#' BrainSpace
#' 
#' This class represents the geometry of a brain image
#' @rdname BrainSpace-class
#' @slot Dim the grid dimensions of the image
#' @slot origin the coordinates of the spatial origin
#' @slot spacing the dimensions (in mm) of the grid units (voxels)
#' @slot axes the set of named spatial axes
#' @slot trans an affine transformation matrix that moves from grid -> real world coordinates
#' @slot inverseTrans an inverse matrix that moves from real world -> grid coordinates
#' @export
setClass("BrainSpace",
		representation=
				representation(Dim = "integer", origin = "numeric", spacing = "numeric",
                   axes="AxisSet", trans="matrix", inverseTrans="matrix"),
 
    validity = function(object) {
      Dim <- object@Dim
      if (length(Dim) < length(object@spacing)) {
        return("Dim slot must be of same length as spacing slot")
      }
      if (length(Dim) < length(object@origin)) {
        return("Dim slot must be of same length as origin slot")
      }
      if (length(Dim) < ndim(object@axes)) {
        return("Dim slot must be of same length as number of axes in AxisSet")
      }
      
      if (any(Dim) < 0) {
        return("Dim slot must contain non-negative values")
      }
    })
         


#' BrainData
#' 
#' Base class for brain image data
#' 
#' 
#' @rdname BrainData-class
#' @slot source an instance of class \code{\linkS4class{BaseSource}} to store the source of the data
#' @slot space an instance of class \code{\linkS4class{BrainSpace}} to represent the geometry of the data space
#' @export
setClass("BrainData",
    representation=
			representation(source="BaseSource", 
					       space="BrainSpace"), 
	contains=c("VIRTUAL"))

#' BrainSlice
#' 
#' Two-dimensional brain image
#' @rdname BrainSlice-class
#' @export
setClass("BrainSlice",       
	    contains=c("BrainData", "array"))

#' Base class for image representing 3D volumetric data.
#' @rdname BrainVolume-class
#' @export
setClass("BrainVolume", 	
	    contains=c("BrainData"))



#' DenseBrainVolume
#' 
#' Three-dimensional brain image, backed by an \code{array}	   
#' @rdname DenseBrainVolume-class
#' @export 
setClass("DenseBrainVolume", 	
		contains=c("BrainVolume", "array"))


#' SparseBrainVolume
#' 
#' Three-dimensional brain image, backed by a \code{sparseVector} for \code{Matrix} package
#' @slot data a \code{sparseVector} instance
#' @rdname SparseBrainVolume-class
#' @export 
setClass("SparseBrainVolume",   
         representation=representation(data="sparseVector"),
         contains=c("BrainVolume"))


#' LogicalBrainVolume
#' 
#' Three-dimensional brain image where all values are either TRUE or FALSE	   
#' @rdname LogicalBrainVolume-class
#' @export  
setClass("LogicalBrainVolume", 	
		contains=c("DenseBrainVolume"))

#' ClusteredBrainVolume
#' 
#' Three-dimensional brain image that is divided into N disjoint partitions    
#' @rdname ClusteredBrainVolume-class
#' @export
setClass("ClusteredBrainVolume",   
         representation=representation(mask="LogicalBrainVolume", clusters="integer", centers="matrix", labelMap="list", clusterMap="hash"),
         contains=c("BrainVolume"))


#' IndexLookupVolume
#' 
#' Three-dimensional brain image that can be used as a map between 1D grid indices and a table of values
#' Currently used in the \code{\linkS4class{SparseBrainVector}} class. 
#' @rdname IndexLookupVolume-class
#' @export
setClass("IndexLookupVolume", 
		representation=
				representation(space="BrainSpace", indices="integer", map="integer"),
		contains=c("BrainVolume"))



#' BrainVector
#' 
#' Four-dimensional brain image	   
#' 
#' @rdname BrainVector-class
#' @export
setClass("BrainVector", 
		contains=c("BrainData"))

#' DenseBrainVector
#' 
#' Four-dimensional brain image, backed by an array   
#' @name DenseBrainVector-class
#' @export 
setClass("DenseBrainVector", 
		contains=c("BrainVector", "array"))

# DenseMMapBrainVector
# 
# Four-dimensional brain image, backed by a memory-mapped file
# @exportClass DenseMMapBrainVector  
# setClass("DenseMMapBrainVector",
# 		representation=representation(data="mmap"),
#		contains=c("BrainVector"))


#' SparseBrainVector
#' 
#' a sparse four-dimensional brain image, backed by a \code{matrix}, where each column represents 
#' a vector spanning the fourth dimension (e.g. time)
#' @rdname SparseBrainVector-class
#' @slot mask the mask defining the sparse domain
#' @slot data the matrix of series, where rows span across voxel space and columns span the fourth dimensions
#' @slot map instance of class \code{\linkS4class{IndexLookupVolume}} is used to map between spatial and index/row coordinates
#' @export
setClass("SparseBrainVector", 
		representation=representation(mask="LogicalBrainVolume",data="matrix", map="IndexLookupVolume"),
		contains=c("BrainVector")) 

#' SparseBrainVectorSource
#' 
#' A class that is used to produce a \code{\linkS4class{SparseBrainVector}} instance
#' @rdname SparseBrainVectorSource-class
#' @slot mask the subset of voxels that will be stored in memory
#' @export
setClass("SparseBrainVectorSource", representation=
				representation(mask="LogicalBrainVolume"),
		contains=c("BrainVectorSource"))



#setClass("TiledBrainVector", 		 
#		representation=representation(cache="list", filename="character", indexList="list", mask="BrainVolume",capacity="numeric"),		
#	    contains=c("BrainVector"))

#' ROIVolume
#' 
#' A class that representing a volumetric region of interest (ROI).
#' @rdname ROIVolume-class
#' @slot data the \code{numeric} data stored in the ROI
#' @slot coords the voxel coordinates of the ROI
#' @exportClass ROIVolume
setClass("ROIVolume", 
		representation=representation(data="numeric", coords="matrix"), contains=c("BrainData"),
         validity = function(object) {
           if (ncol(object@coords) != 3) {
             stop("coords slot must be a matrix with 3 columns")
           }
         })

#' Kernel
#' 
#' A class representing an image kernel
#' 
#' @rdname Kernel-class
#' @slot width the width in voxels of the kernel
#' @slot weights the kernel weights
#' @slot voxels the relative voxel coordinates of the kernel
#' @slot coords the relative real coordinates of the kernel
#' @export
setClass("Kernel", 
         representation=representation(width="numeric", weights="numeric", voxels="matrix", coords="matrix"))
         

## TODO add a LazyBrainBucket class

#' BrainSurface
#' 
#' a three-dimensional surface consisting of a set of triangle vertices with one value per vertex.
#' @rdname BrainSurface-class
#' @slot source the data source for the surface
#' @slot mesh the underlying \code{mesh3d} object 
#' @slot data the vector of data value at each vertex of the mesh
#' @export
setClass("BrainSurface", 
         representation=representation(source="BrainSource", mesh="mesh3d", data="numeric"))

#' BrainSurfaceVector
#' 
#' a three-dimensional surface consisting of a set of triangle vertices with multiple values per vertex.
#' 
#' @rdname BrainSurfaceVector-class
#' @slot source the data source for the surface
#' @slot mesh the underlying \code{mesh3d} object 
#' @slot mat a matrix of values where each column contains a vector of values over the surface nodes.
#' @export
setClass("BrainSurfaceVector", 
         representation=representation(source="BrainSource", mesh="mesh3d", mat="numeric"))

#' BrainBucket
#' 
#' a four-dimensional image that conists of a sequence of labeled image volumes backed by a list
#' @rdname BrainBucket-class
#' @slot source the data source for the bucket volumes
#' @slot labels the names of the sub-volumes contained in the bucket
#' @slot data a list of \code{\linkS4class{BrainVolume}} instances with names corresponding to volume labels
#' @export
setClass("BrainBucket", 
		representation=representation(source="BrainSource", labels="character", data="list"),
		validity = function(object) {
		  if (any(sapply(object@data, function(obj) !is(obj, "BrainVolume")))) {
        stop("all elements of data list must be of type `BrainVolume`")
		  } else {
        TRUE
		  }
		},
		contains=c("BrainVector"))
				

#setClassUnion(name="index", members =  c("numeric", "logical", "character"))

#' Layer
#' 
#' A class used for displaying 2D images with color maps
#' 
#' @rdname Layer-class
#' @slot vol the \code{BrainVolume} that provides the data for the layer.
#' @slot colorMap a character vector of colors in hexadecimal rgb format. 
#'       Can be generated by calls to \code{rainbow}, \code{heat.colors}, \code{topo.colors}, \code{terrain.colors} or similar functions.
#' @slot thresh cut-off value above which vlaues will be made transparent.
#' @slot axis the axis index of perpendicular to the xy plane (option: 1,2,3; default is 3)
#' @slot zero.col the color pixels with intensity of zero. This value overrides the color from the slot \code{colorMap}
#' @slot alpha the transparency of the layer
#' @export
setClass("Layer",
         representation=representation(vol="BrainVolume", colorMap="vector", thresh="numeric", axis="numeric", zero.col="character", alpha="numeric"))

setClass("Overlay",
                  representation(layers="list"))
                  

setClass("RenderedSlice",
         representation=representation(slice="BrainSlice",
                                       width="numeric",
                                       height="numeric",
                                       raster="rastergrob"))

setClass("RenderedSliceStack",
         representation=representation(slices="list",
                                       width="numeric",
                                       height="numeric",
                                       grob="gList"))
