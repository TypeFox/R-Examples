# Author: Shaun Walbridge
# prevent directory colisions on multiuser machines by generating a unique dir

.tmppath <- function() {

# default temp files now inside the session tempdir. Therefore no need to add username etc.
   return( file.path(tempdir(), 'raster', '/') )

   
# when Sys.info is NULL, use this default
	extension <- 'user'
	s <- Sys.info()
	if (!is.null(s)) {
		# get userid from system, to generate temporary directory name
		user <- s[["user"]]
		if (!is.null(user)) {
			extension <- user
		}
	}
#	d <- paste(dirname(tempdir()), '/R_raster_tmp/', extension, '/', sep="")
    d <- paste(dirname(tempdir()), '/R_raster_', extension, '/', sep="")	
	return(d)
}

