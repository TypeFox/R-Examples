# Functions to identify and read a variety of bitmap image formats.
###############################################################################

#' Identify the type of an image using the magic value at the start of the file 
#'
#' Currently works for png, jpeg and BMP images.
#' Will seek to start of file if passed a connection.
#' For details of magic values for files, see e.g. 
#' http://en.wikipedia.org/wiki/Magic_number_(programming)#Magic_numbers_in_files
#' @param source Path to file or connection
#' @param Verbose Whether to write a message to console on failure (Default F)
#' @return character value corresponding to standard file extension of 
#'   image format (i.e. jpg, png, bmp) or NA_character_ on failure.
#' @export
#' @examples
#' jpegfile=system.file("img", "Rlogo.jpg", package="jpeg")
#' image_type(jpegfile)
#' jpeg_pretending_to_be_png=tempfile(fileext = '.png')
#' file.copy(jpegfile, jpeg_pretending_to_be_png)
#' image_type(jpeg_pretending_to_be_png)
#' unlink(jpeg_pretending_to_be_png)
image_type<-function(source,Verbose=FALSE){
  if (inherits(source, "connection")) 
    seek(source, 0)
  magic = readBin(source, what = 0L, n = 8, size = 1L, signed = FALSE)
  if(isTRUE(all.equal(magic[1:2], c(66, 77))))
    return('bmp')
  else if(isTRUE(all.equal(magic[1:8], 
          c(0x89,0x50,0x4E,0x47,0x0D, 0x0A, 0x1A, 0x0A))) )
    return('png')
  else if(isTRUE(all.equal(magic[1:2], c(0xFF, 0xD8))))
    return('jpg')
  # otherwise we failed to identify the file
  if(Verbose) warning("Failed to identify image type of: ",source,
        ' with magic: ',format.hexmode(as.raw(magic)))
  return(NA_character_)
}

#' Read in a bitmap image in JPEG, PNG or BMP format
#'
#' By default uses magic byte to identify file 
#'   (rather than the file extension)
#' Currently uses readers in bmp, jpeg and png packages.
#' @param f Path to image file
#' @param channel Integer identifying channel to return for an RGB image
#' @param IdentifyByExtension Identify by file extension only (Default FALSE)
#' @param ... Additional parameters passed to underlying image readers
#' @return return value
#' @importFrom png readPNG
#' @importFrom jpeg readJPEG
#' @importFrom bmp read.bmp
#' @export
#' @seealso \code{\link[jpeg]{readJPEG},\link[png]{readPNG},\link[bmp]{read.bmp}}
#' @examples
#' img1=read.bitmap(system.file("img", "Rlogo.jpg", package="jpeg"))
#' str(img1)
#' img2 <- read.bitmap(system.file("img", "Rlogo.png", package="png"))
#' # nb the PNG image has an alpha channel
#' str(img2)
read.bitmap<-function(f,channel,IdentifyByExtension=FALSE,...){
  
  if(!file.exists(f)) stop("File: ",f," does not exist.")
  
  if(IdentifyByExtension) 
    ext=tolower(sub(".*\\.([^.]+)$","\\1",f))
  else
    ext=image_type(f)
  
  readfun=switch(ext,png=readPNG,jpeg=readJPEG,jpg=readJPEG,bmp=read.bmp,
      stop("File f: ",f," does not appear to be a PNG, BMP or JPEG"))
  im=readfun(f,...)
  
  if(!missing(channel) && length(dim(im))==3) im=im[,,channel]
  im
}
