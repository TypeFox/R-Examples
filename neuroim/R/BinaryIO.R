
#' @include AllClass.R
roxygen()
#' @include AllGeneric.R
roxygen()

#' ColumnReader
#' 
#' Constructor for  \code{\linkS4class{ColumnReader}} class
#' 
#' @param nrow the number of rows
#' @param ncol the number of columns
#' @param reader a function that takes a set of column indices and returns a \code{matrix}
#' @rdname ColumnReader
#' @export 
ColumnReader <- function(nrow, ncol, reader) {
  stopifnot(is.function(reader))
  new("ColumnReader", nrow=as.integer(nrow), ncol=as.integer(ncol), reader)
}



#' BinaryReader
#' 
#' Constructor for  \code{\linkS4class{BinaryReader}} class
#' 
#' @param input file name to read from or else a \code{connection} object
#' @param byteOffset the number of bytes to skip at the start of input
#' @param dataType R data type of binary elements
#' @param bytesPerElement number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @param endian endianness of binary input connection
#' @rdname BinaryReader
#' @export 
BinaryReader <- function(input, byteOffset, dataType, bytesPerElement, endian=.Platform$endian) {
	if (is.character(input)) { 
		new("BinaryReader", input=file(input, open="rb"), byteOffset=as.integer(byteOffset), dataType=dataType, bytesPerElement=as.integer(bytesPerElement), endian=endian)
	} else {
		stopifnot(inherits(input, "connection")) 
		new("BinaryReader", input=input, byteOffset=as.integer(byteOffset), dataType=dataType, bytesPerElement=as.integer(bytesPerElement), endian=endian)		
	}

}

#' BinaryWriter
#' 
#' Constructor for  \code{\linkS4class{BinaryWriter}} class
#' @param output file name to write to or else a \code{connection} object
#' @param byteOffset the number of bytes to skip at the start of output
#' @param dataType R data type of binary elements
#' @param bytesPerElement number of bytes in each data element (e.g. 4 or 8 for floating point numbers)
#' @param endian endianness of binary output connection
#' @rdname BinaryWriter-class
#' @export 
BinaryWriter <- function(output, byteOffset, dataType, bytesPerElement, endian=.Platform$endian) {
	if (is.character(output)) { 
		new("BinaryWriter", output=file(output, open="wb"), byteOffset=as.integer(byteOffset), dataType=dataType, bytesPerElement=as.integer(bytesPerElement), endian=endian)
	} else {
		stopifnot(inherits(output, "connection")) 
		new("BinaryWriter", output=output, byteOffset=as.integer(byteOffset), dataType=dataType, bytesPerElement=as.integer(bytesPerElement), endian=endian)		
	}
	
}

## code duplication, fix me. introduce "BinaryConnection superclass
setMethod(f="initialize", signature=signature("BinaryReader"), 
		def=function(.Object, input, byteOffset, dataType, bytesPerElement, endian) {
			.Object@input <- input
			.Object@byteOffset <- byteOffset
			.Object@dataType <- dataType
			.Object@bytesPerElement <- bytesPerElement
			.Object@endian <- endian
			
			## must be seekable connection, should enforce this
			## 
			
			if (attr(.Object@input, "class")[[1]] != "gzfile") {
				seek(.Object@input, where=.Object@byteOffset, origin="start")	
			} else {		
				n <- as.integer(.Object@byteOffset/.Object@bytesPerElement)
				readBin(.Object@input, what=.Object@dataType, size=.Object@bytesPerElement, endian=.Object@endian, n=n)
			}
			
			.Object
			
		}) 

## code duplication, fix me
setMethod(f="initialize", signature=signature("BinaryWriter"), 
		def=function(.Object, output, byteOffset, dataType, bytesPerElement, endian) {
			.Object@output <- output
			.Object@byteOffset <- byteOffset
			.Object@dataType <- dataType
			.Object@bytesPerElement <- bytesPerElement
			.Object@endian <- endian
			
			## must be seekable connection, should enforce this
			## 
			#seek(.Object@output, where=.Object@byteOffset, origin="start")		
			.Object
		})  

#' readElements
#' 
#' @export
#' @rdname readElements-methods
setMethod(f="readElements", signature=signature(x= "BinaryReader", numElements="numeric"),
		def=function(x, numElements) {
			readBin(x@input, what=x@dataType, size=x@bytesPerElement, n=numElements, endian=x@endian)		
		})

#' writeElements
#' 
#' @rdname writeElements-methods
setMethod(f="writeElements", signature=signature(x= "BinaryWriter", els="numeric"),
		def=function(x, els) {
			if (.getRStorage(x@dataType) == "integer") {
				writeBin(as.integer(els), x@output, size=x@bytesPerElement, endian=x@endian)		
			} else if (.getRStorage(x@dataType) == "double") {
				writeBin(as.double(els), x@output, size=x@bytesPerElement, endian=x@endian)	
			} else {
				stop("failed to handle data type: ", x@dataType)
			}				
			
		})

#' close
#'  
#' @rdname close-methods
#' @export
#' @param con the object to close
setMethod(f="close", signature=signature(con= "BinaryReader"),
		def=function(con) {
			base::close(con@input)				
		})


#' @rdname close-methods
setMethod(f="close", signature=signature(con= "BinaryWriter"),
		def=function(con) {
			base::close(con@output)				
		})
