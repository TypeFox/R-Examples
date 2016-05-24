#' Extract the CRC (32 bit hash) of a gzip file
#' 
#' @description  Reads the crc from a gzip file, assuming it is the last 4 bytes
#'   of the file. First checks for a valid gzip magic number at the start of the
#'   file.
#' @details CRC32 is not a strong hash like SHA1 or even MD5, but it does 
#'   provide a basic hash of the \bold{uncompressed contents} of the gzip file.
#'   NB CRCs are stored in little endian byte order regardless of platform.
#' @param f Path to a gzip file
#' @return hexadecimal formatted
#' @export
#' @examples
#' rdsfile=system.file('help/aliases.rds')
#' gzip.crc(rdsfile)
gzip.crc<-function(f){
  con=file(f,open='rb')
  on.exit(close(con),add=TRUE)
  magic=readBin(con,what='raw',n=2)
  if(magic[1]!=0x1f || magic[2]!=0x8b) {
    warning("This is not a gzip file")
    return(NA)
  }
  seek(con,-8,origin='end')
  # nb gzip always writes little endian
  crc=readBin(con,integer(),size=4, endian = 'little')
  format(as.hexmode(crc),width=8)
}

#' Check if a file is a gzip file
#' 
#' @param f Path to file to test
#' @return logical indicating whether \code{f} is in gzip format (or \code{NA}
#'   if the file cannot be accessed)
#' @export
#' @examples
#' notgzipfile=tempfile()
#' writeLines('not a gzip', notgzipfile)
#' is.gzip(notgzipfile)
#' con=gzfile(gzipfile<-tempfile(),open='wt')
#' writeLines('This one is gzipped', con)
#' close(con)
#' is.gzip(gzipfile)
#' unlink(c(notgzipfile,gzipfile))
is.gzip<-function(f) {
  if(file.access(f, mode=4)<0) return(NA)
  x=file(f)
  on.exit(close(x))
  isTRUE(summary(x)$class=='gzfile')
}
