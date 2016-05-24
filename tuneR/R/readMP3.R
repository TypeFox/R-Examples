## decode.R - decode MP3 files
##
## Author: Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>

readMP3 <- function(filename){
    if(!is.character(filename))
        stop("'filename' must be of type character.")
    if(length(filename) != 1)
        stop("Please specify exactly one 'filename'.")
    if(!file.exists(filename))
        stop("File '", filename, "' does not exist.")
    if(file.access(filename, 4))
        stop("No read permission for file ", filename)
    con <- file(filename, "rb")
    on.exit(close(con)) # be careful ...

    data <- readBin(con, raw(), n = file.info(filename)$size)
    .Call("do_read_mp3", data, package = "tuneR")
}
