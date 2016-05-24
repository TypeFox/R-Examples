writeFITSim <-
function (X, file = "R.fits", type = "double", bscale = 1,
    bzero = 0, c1 = NA, c2 = NA, crpixn = NA, crvaln = NA, cdeltn = NA,
    ctypen = NA, cunitn = NA, axDat = NA, header = '')
{
### Function writes FITS image after primary header
###
### Takes:
  ## Multi-dimensional array: X
  ## Output file name: file
  ## Type to write (single or double precision)
  ## Overall scaling and shifting variables: bscale, bzero
  ## Text comment variables: c1 and c2
  ## Axis parameters with usual FITS standard meanings
### Returns:
  ## Writes FITS file to disk
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001)
###
### A. Harris, Univ. MD Astronomy, 12/30/2012
###   (adapted from version of 4/22/08, 4/11/12, retaining items for
###    backwards compatibility)

    ## Open file
    zz <- file(file, "wb")

    ## Generate and write a primary header
    hdr0 <- makeFITSimHdr(X, primaryhdu = TRUE, type = type,
                          c1 = c1, c2 = c2, bscale = bscale, bzero = bzero,
                          crpixn = crpixn, crvaln = crvaln,
                          cdeltn = cdeltn, ctypen = ctypen, cunitn = cunitn,
                          axDat = axDat, header = header)
    writeChar(hdr0, zz, eos = NULL)

    ## Establish size to write, write data
        ## Choose integer, float; byte single, double precision based on data
    type <- tolower(substr(type, 1, 1)) # double, single, byte
    if (is.integer(X)) {
        switch(type, b = {
            size <- 1
        }, s = {
            size <- 2
        }, d = {
            size <- 4
        }, stop("Unrecognized data type: not single, double, or byte"))
    }
    else {
        switch(type, s = {
            size <- 4
        }, d = {
            size <- 8
        }, stop("Unrecognized data type: not single or double"))
    }
    writeBin(as.vector(X), zz, size = size, endian = "big")

    ## Pad rest of record with zeros and close
    pad <- raw(2880 - (length(as.vector(X)) * size)%%2880)
    writeBin(pad, zz, endian = "big")
    close(zz)
}


