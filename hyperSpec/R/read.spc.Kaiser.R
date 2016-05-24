##' Import functions for Kaiser Optical Systems .spc files
##' 
##' \code{read.spc.Kaiser} imports sets of .spc files written by Kaiser Optical Systems' Hologram
##' software.  It may also serve as an example how to write wrapper functions for \code{read.spc} to
##' conveniently import specialized sets of .spc files.
##' 
##' @title read Kaiser .spc files
##' @export
##' @rdname read-spc-Kaiser
##' @param files If \code{glob = TRUE}, \code{filename} can contain wildcards.
##'   Thus all files matching the name pattern in \code{filename} can be
##'   specified.
##' @param glob If \code{TRUE} the filename is interpreted as a wildcard
##'   containing file name pattern and expanded to all matching file names.
##' @param keys.log2data,... All further arguments are handed over directly to \code{\link{read.spc}}.
##' @return hyperSpec 
##' @examples
##' ## for examples, please see `vignette ("fileio", package = "hyperSpec")`.

read.spc.Kaiser <- function (files, ..., glob = TRUE) {
    
	if (glob)
		files <- Sys.glob (files)

   if (length (files) == 0){
     warning ("No files found.")
     return (new ("hyperSpec"))
   }
	
	f <- files [1]
	
	spc <- read.spc (f, no.object = TRUE, ...)
	
	data <- spc$data [rep (1L, length (files)),, drop = FALSE]
	
	spc$spc  <- spc$spc  [rep (1L, length (files)), , drop = FALSE]
	
	for (f in seq_along (files)){
		tmp <- read.spc (files [f], no.object = TRUE, ...)
		
		data [f, ] <- tmp$data 
		spc$spc  [f, ] <- tmp$spc
	}

	data$file <- files
   
	new ("hyperSpec", wavelength = spc$wavelength, spc = spc$spc, data = data, 
			labels = tmp$label)
}

##' \code{read.spc.KaiserMap} is a wrapper for \code{read.spc.Kaiser} with predefined \code{log2data}
##' to fetch the stage position for each file.
##' @rdname read-spc-Kaiser
##' @export
read.spc.KaiserMap <- function (files, keys.log2data = NULL, ...) {
  keys.log2data <- c ('Stage_X_Position','Stage_Y_Position','Stage_Z_Position', keys.log2data)

  spc <- read.spc.Kaiser (files, keys.log2data = keys.log2data, ...)
  
  spc@data <- spc@data [, ! colnames (spc@data) %in% c ("z", "z.end"), drop = FALSE]

  colnames (spc@data) <- gsub ("Stage_(.)_Position", "\\L\\1", colnames (spc@data), perl = TRUE)
  for (cln in c ("x", "y", "z"))
      spc@data [[cln]] <- as.numeric (spc@data [[cln]])

  spc@label$x <- expression (`/` (x, micro * m))
  spc@label$y <- expression (`/` (y, micro * m))
  spc@label$z <- expression (`/` (z, micro * m))
  spc@label$z.end <- NULL
  
  spc
}

##' \code{read.spc.KaiserLowHigh} is a wrapper for \code{read.spc.Kaiser} for raw data that is saved
##' in separate files for low and high wavenumber range.  The wavelength axis holds the pixel
##' numbers, which repeat for low and high wavenumber ranges. 
##' 
##' @rdname read-spc-Kaiser
##' @param type what kind of measurement was done? If \code{"map"}, \code{read.spc.KaiserMap} is used
##' instead of \code{read.spc.Kaiser}.
##' @export
read.spc.KaiserLowHigh <- function (files = stop ("file names needed"),
                                    type = c ("single", "map"),
                                    ..., glob = TRUE) {

	if (glob)
       files <- Sys.glob (files)
   
   files <- matrix (files, nrow = 2)

   type <- match.arg (type)
   switch (type,
           single = cbind (read.spc.Kaiser    (files [1,], ..., glob = FALSE),
                           read.spc.Kaiser    (files [2,], ..., glob = FALSE)),
           map    = cbind (read.spc.KaiserMap (files [1,], ..., glob = FALSE),
                           read.spc.KaiserMap (files [2,], ..., glob = FALSE))
           )

}
