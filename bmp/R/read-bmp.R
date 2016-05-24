# Functions to read Windows BMP file format
# 
# Author: jefferis
###############################################################################

# Header format
# Data	Description
# WORD Type;	File type. Set to "BM".
# DWORD Size;	Size in BYTES of the file.
# DWORD Reserved;	Reserved. Set to zero.
# DWORD Offset;	Offset to the data.
# DWORD headerSize;	Size of rest of header. Set to 40.
# DWORD Width;	Width of bitmap in pixels.
# DWORD Height;	Height of bitmap in pixels.
# WORD Planes;	Number of Planes. Set to 1.
# WORD BitsPerPixel;	Number of bits per pixel.
# DWORD Compression;	Compression. Usually set to 0.
# DWORD SizeImage;	Size in bytes of the bitmap.
# DWORD XPixelsPerMeter;	Horizontal pixels per meter.
# DWORD YPixelsPerMeter;	Vertical pixels per meter.
# DWORD ColorsUsed;	Number of colors used.
# DWORD ColorsImportant;	Number of "important" colors.

#' Open windows BMP format image files
#' 
#' Limited to 8 bit greyscale images and 24 bit RGB images.
#' @param f File to open
#' @param Verbose Give verbose warnings (default FALSE)
#' @return array of dims height x width x channels 
#' @author jefferis
#' @export
#' @examples
#' \dontrun{
#' library(pixmap)
#' r=read.bmp('myrgbimage.bmp')
#' pr=pixmapRGB(r)
#' r=read.bmp('mygreyimage.bmp')
#' pr=pixmapGrey(r)
#' plot(pr)
#' } 
read.bmp<-function(f,Verbose=FALSE){
	con=file(f,open='rb')
  on.exit(close(con))
  if(!is.bmp(con))
    stop(basename(f)," is not a valid BMP file")
  # file header
  h=list()
  h$filesize=ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  readBin(con,what=1L,size=4,endian='little')
  h$offset=ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  # image header
  h$header_sz=ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  h$width=ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  h$height=ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  h$nplanes=readBin(con,what=1L,size=2,endian='little',signed = FALSE)
  h$depth = readBin(con,what=1L,size=2,endian='little',signed = FALSE)
  h$compress_type = ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  h$bmp_bytesz = ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  h$hres = ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  h$vres = ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  h$ncolors = ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  h$nimpcolors = ConvertIntToUInt(readBin(con,what=1L,size=4,endian='little'))
  
  if(h$compress_type!=0)
     stop("Do not know how to handle compressed BMP")
  if(!h$depth %in% c(8,24,32))
    stop("Do not know how to handle bit depth: ",h$depth)
  
  bytes_pixel=h$depth / 8
  row_width_bytes=h$width * bytes_pixel
  # bmp rows are written to lenth of nearest 4 bytes
  rounded_row_width_bytes = ceiling(row_width_bytes/4)*4
  bytes_to_trim = row_width_bytes %% 4
  bytes_to_read=rounded_row_width_bytes * h$height
  if(h$bmp_bytesz==0) {
    if(Verbose) warning("invalid byte size information for image")
  }
  else if(h$bmp_bytesz != bytes_to_read)
    stop("mismatch between predicted and actual number of bytes in image")
  seek(con,h$offset)
  if(h$depth>8){
    nchans=h$depth/8
    d=readBin(con,what=1L,size=1,n=bytes_to_read,endian='little',signed=FALSE)
    dim(d)=c(rounded_row_width_bytes,h$height)
    # trim off padding bytes
    d=d[1:row_width_bytes,]
    # dims on read are colours, columns, rows
    dim(d)=c(nchans,h$width,h$height)
    # swap row and colour dimensions
    d=aperm(d,c(3,2,1))
    # flip rows (i.e. Y) and colours (since little endian => BGR)
    d=d[h$height:1,,nchans:1]
  }
  else {
    d=readBin(con,what=1L,size=bytes_pixel,n=bytes_to_read/bytes_pixel,
        endian='little',signed=FALSE)
    dim(d)=c(rounded_row_width_bytes / bytes_pixel,h$height)
    # trim off excess bytes and transpose
    d=t(d[1:(row_width_bytes / bytes_pixel),])
    # flip up/down (bitmaps start in bottom left)
    d = d[nrow(d):1,]
  }
  attr(d,'header')=h
  
  d
}

#' Returns TRUE if file is a Windows BMP image
#' 
#' NB this just checks the magic 'BM' in the first two bytes of the file
#' @param source file or connection
#' @return TRUE or FALSE
#' @author jefferis
#' @export
is.bmp<-function(source) {
  if(inherits(source,'connection'))
    seek(source,0)
	magic=readBin(source,what=0L,n=2,size=1L)
  isTRUE(all.equal(magic,c(66,77)))
}

#' Fix a 32 bit unsigned integer that has been read as signed
#' 
#' This is really just to fix a limitation of readBin/R's 32 bit signed ints  
#' @param x Number to be fixed 
#' @param adjustment number to be added to convert to uint32 (2^32 by default)
#' @return numeric value of uint32
#' @author jefferis
#' @seealso \code{\link{readBin}}
ConvertIntToUInt<-function(x,adjustment=2^32){
	x=as.numeric(x)
	signs=sign(x)
	x[signs<0]=x[signs<0]+adjustment
	x
}
