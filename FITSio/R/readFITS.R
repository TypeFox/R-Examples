readFITS <-
function (file = "R.fits", hdu = 1, maxLines = 5000,
          fixHdr = c('none', 'remove', 'substitute'), phdu = 1)
{
### Simple reader for FITS bintable, array, and image files
### Also serves as template for mulitple header and data unit reads
###
### Takes:
  ## FITS file name: file
  ## Read nth header and data unit from file: hdu
  ## Set phdu=0 in call if NAXIS!=0 but header is in secondary header unit
### Returns:
  ## List with data, parameters, and header vector
### Requires/Used by:
  ## Requires readFITSbintable.r
  ## Requires readFITSarray.r
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001)
###
### A. Harris, Univ. MD Astronomy, 4/22/08
  ## Added multiple image reads 9/22/10 AH
  ## Updated for full header and new header parsing 12/31/12 AH
###
    ## Open file, read primary header unit
    zz <- file(file, "rb")
    header <- readFITSheader(zz, maxLines = maxLines, fixHdr = fixHdr[1])
    hdr <- parseHdr(header)

    ## Determine number of array dimensions, select appropriate extension
    tmp <- hdr[which(hdr == "NAXIS") + 1]
    if (tmp == "")
        tmp <- "0"
    if (as.numeric(tmp) > 0) {
        tmp <- as.numeric(hdr[which(hdr == "NAXIS1") + 1])
        if (tmp == 0)
            phdu <- 0
    }
    if ((as.numeric(tmp) * phdu) > 0) {  # Data are an array off primary header
        D <- readFITSarray(zz, hdr)
        if (hdu > 1) {
            for (i in 2:hdu) {    # Brute-force read to hdu-th header-data unit
                header <- readFITSheader(zz, maxLines = maxLines, fixHdr = fixHdr[1])
                hdr <- parseHdr(header)
                switch(tolower(hdr[which(hdr == "XTENSION") + 1]),
                       bintable = {
                           D <- readFITSbintable(zz, hdr)
                       }, image = {
                           D <- readFITSarray(zz, hdr)
                       }, stop("Current version supports only bintable and image"))
            }
        }
        close(zz)
        D$header <- header
        return(D)
    }
    else {                    # Data are off an extension header
        for (i in 1:hdu) {    # Brute-force read to hdu-th header-data unit
            header <- readFITSheader(zz, maxLines = maxLines, fixHdr = fixHdr[1])
            hdr <- parseHdr(header)
            switch(tolower(hdr[which(hdr == "XTENSION") + 1]),
                bintable = {
                  D <- readFITSbintable(zz, hdr)
                }, image = {
                  D <- readFITSarray(zz, hdr)
                }, stop("Current version supports only bintable and image"))
        }
        close(zz)
        D$header <- header
        return(D)
    }
}

