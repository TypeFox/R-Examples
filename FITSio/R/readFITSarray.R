readFITSarray <-
function (zz, hdr)
{
### Reader for FITS multidimentsional arrays, including 2D images
###
### Takes:
  ## File handle: zz
  ## Parsed header vector: hdr
### Returns:
  ## Data array: imDat
  ## Axis data frame: axDat
  ## Header vector: hdr
### Requires/Used by:
  ## Requires readFITSheader.r
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001)
###
### A. Harris, Univ. MD Astronomy, 4/17/08
### Changed defaults for missing CRPIX, CRVAL, and CDELT from NA to 1 (L65-71)
### AH 7/11/09
### Fixed apparent problem in padding calculation and read, 9/28/10 AH
### Updated for new header handling, 12/30/12 AH

    ## Parse header if full header is supplied instead of parsed version
    if (nchar(hdr[1])==80) hdr <- parseHdr(hdr)

    ## Determine number of array dimensions
    naxis <- as.numeric(hdr[which(hdr == "NAXIS") + 1])
    ## Find the right data type
    switch(hdr[which(hdr == "BITPIX") + 1], "-64" = {
        bsize <- 8
        btype <- numeric()
        bsign <- TRUE
    }, "-32" = {
        bsize <- 4
        btype <- numeric()
        bsign <- TRUE
    }, "32" = {
        bsize <- 4
        btype <- integer()
        bsign <- TRUE
    }, "16" = {
        bsize <- 2
        btype <- integer()
        bsign <- TRUE
    }, "8" = {
        bsize <- 1
        btype <- integer()
        bsign <- FALSE
    }, stop("Unknown BITPIX request in readFITSarray"))
    ## Put array information into vectors
    NAXISn <- integer(naxis)
    CRPIXn <- integer(naxis)
    CRVALn <- numeric(naxis)
    CDELTn <- numeric(naxis)
    CTYPEn <- character(naxis)
    CUNITn <- character(naxis)
    PTYPEn <- character(naxis)
    PSCALn <- numeric(naxis)
    PZEROn <- numeric(naxis)
    numwords <- 1
    for (i in 1:naxis) {
        tmp <- as.numeric(hdr[which(hdr == paste("NAXIS", i,
            sep = "")) + 1])
        numwords <- numwords * tmp
        NAXISn[i] <- tmp
        tmp <- hdr[which(hdr == paste("CRPIX", i, sep = "")) +
            1]
        CRPIXn[i] <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
        tmp <- hdr[which(hdr == paste("CRVAL", i, sep = "")) +
            1]
        CRVALn[i] <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
        tmp <- hdr[which(hdr == paste("CDELT", i, sep = "")) +
            1]
        CDELTn[i] <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
        tmp <- hdr[which(hdr == paste("CTYPE", i, sep = "")) +
            1]
        CTYPEn[i] <- ifelse(length(tmp) != 1, "", tmp)
        tmp <- hdr[which(hdr == paste("CUNIT", i, sep = "")) +
            1]
        CUNITn[i] <- ifelse(length(tmp) != 1, "", tmp)
    }
    ## Read data into array.
    D <- array(readBin(zz, what = btype, n = numwords, size = bsize,
                       signed = bsign, endian = "big"), dim = NAXISn)
    ## Finish reading block
    nbyte <- numwords * bsize
    nbyte <- ifelse(nbyte%%2880 == 0, 0, 2880 - nbyte%%2880)
    tmp <- readBin(zz, what = 'raw', n = nbyte)
    ## Scale image to physical units if needed
    tmp <- hdr[which(hdr == "BUNIT") + 1]
    BUNIT <- ifelse(length(tmp) != 1, "", tmp)
    tmp <- hdr[which(hdr == "BSCALE") + 1]
    BSCALE <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
    tmp <- hdr[which(hdr == "BZERO") + 1]
    BZERO <- ifelse(length(tmp) != 1, 0, as.numeric(tmp))
    if (BSCALE != 1 || BZERO != 0)
        D <- D * BSCALE + BZERO
    ## Make data frame with axis data
    axDat <- data.frame(crpix=CRPIXn, crval=CRVALn, cdelt=CDELTn,
                        len=dim(D), ctype=CTYPEn, cunit=CUNITn,
                        stringsAsFactors=FALSE)
    ## Return structure with data and image information
    list(imDat = D, axDat = axDat, hdr = hdr)
}

