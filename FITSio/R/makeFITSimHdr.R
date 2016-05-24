makeFITSimHdr <-
function (X, primaryhdu = TRUE, type = 'double',
          c1 = NA, c2 = NA, bscale = 1, bzero = 0,
          crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA,
          axDat = NA, header = '')
{
### Function assembles FITS primary header for images
###    (multi-dimensional arrays)
###
### Takes:
  ## Most variables have names as defined in FITS reference
  ## Additional comment lines: c1, c2
### Returns:
  ## Header data for writeFITSim.r
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001)
###
### A. Harris, Univ. MD Astronomy, 2/1/13

    ## Number of axes, parse data type (single or double precision)
    naxisn <- dim(X)
    naxis <- length(naxisn)

    # Number of bits per pixel
    type <- tolower(substr(type, 1, 1)) # double, single, byte
    if (is.integer(X)) {
        switch(type, b = {
            bitpix <- 8
        }, s = {
            bitpix <- 16
        }, d = {
            bitpix <- 32
        }, stop("Unrecognized data type: not single, double, or byte"))
    }
    else {
        switch(type, s = {
            bitpix <- -32
        }, d = {
            bitpix <- -64
        }, stop("Unrecognized data type: not single or double"))
    }

    ## Use whole axis header if given
    if (is.data.frame(axDat)) {
        crpixn <- axDat$crpix
        crvaln <- axDat$crval
        cdeltn <- axDat$cdelt
        ctypen <- axDat$ctype
        cunitn <- axDat$cunit
    } else {
    ## Otherwise make defaults
        if (is.na(crpixn[1]))
            crpixn[1:naxis] <- 1
        if (is.na(crvaln[1]))
            crvaln[1:naxis] <- 1
        if (is.na(cdeltn[1]))
            cdeltn[1:naxis] <- 1
        if (is.na(ctypen[1]))
            ctypen[1:naxis] <- ""
        if (is.na(cunitn[1]))
            cunitn[1:naxis] <- ""
    }
    ## Error checking
    if (length(crpixn) != naxis) stop(' *** crpixn vector length error ***')
    if (length(crvaln) != naxis) stop(' *** crvaln vector length error ***')
    if (length(cdeltn) != naxis) stop(' *** cdeltn vector length error ***')
    if (length(ctypen) != naxis) stop(' *** ctypen vector length error ***')
    if (length(cunitn) != naxis) stop(' *** cunitn vector length error ***')

    ##### Make card images for header from values above

    if (primaryhdu) {
        cimages <- sprintf('%-80s',
            "SIMPLE  = T                      / File conforms to FITS standard")
     } else {
        cimages <- sprintf('%-80s',
            "XTENSION= 'IMAGE'                / Image extension")
    }

    cimages <- addKwv('BITPIX', bitpix, 'number of bits per data pixel', cimages)
    cimages <- addKwv('NAXIS', naxis, 'number of data axes', cimages)

    tmp <- character(naxis)
    for (i in 1:naxis) {
        tmp[i] <- newKwv(sprintf('NAXIS%d', i), naxisn[i], 'length of data axis')
    }
    cimages <- c(cimages, tmp)

    cimages <- c(cimages, sprintf('%-80s',
             "EXTEND  = T                      / FITS dataset may contain extensions"))
    cimages <- addComment('  Written by the R language FITSio package', cimages)
    cimages <- addComment('  FITS (Flexible Image Transport System) format is defined in', cimages)
    cimages <- addComment('  Astronomy and Astrophysics, volume 376, page 359 (2001)', cimages)

    tmp <- character(naxis)
    j <- 1
    for (i in 1:naxis) {
        tmp[j] <-   newKwv(sprintf('CRPIX%-3d', i), crpixn[i])
        tmp[j+1] <- newKwv(sprintf('CRVAL%-3d', i), crvaln[i])
        tmp[j+2] <- newKwv(sprintf('CDELT%-3d', i), cdeltn[i])
        tmp[j+3] <- newKwv(sprintf('CTYPE%-3d', i), ctypen[i])
        tmp[j+4] <- newKwv(sprintf('CUNIT%-3d', i), cunitn[i])
        j <- j+5
    }
    cimages <- c(cimages, tmp)

    # Add comments and scaling (to keep backwards compatibility)
    if (!is.na(c1)) cimages <- addComment(c1, cimages)
    if (!is.na(c2)) cimages <- addComment(c2, cimages)
    if (!(isTRUE(all.equal(bscale, 1)) & isTRUE(all.equal(bzero, 0)))) {
        cimages <- addKwv('BSCALE', bscale, 'overall scaling', cimages)
        cimages <- addKwv('BZERO', bzero, 'overall offset', cimages)
    }

    # Prevent duplicated reserved keywords before merging with new card images
    if (length(header) > 0) {
        # reserved keywords that must not be duplicated in cimages
        reserved <- c('^ *SIMPLE ', '^XTENSION=',
                 '^ *BITPIX ', '^ *NAXIS', '^ *EXTEND ',
                 '^ *CRPIX', '^ *CRVAL', '^ *CDELT', '^ *CTYPE', '^ *CUNIT',
                 '^ *BSCALE ', '^ *BZERO ', '^ *END ',
                 '^ *COMMENT   Written by the R language FITSio package',
                 '^ *COMMENT   FITS \\(Flexible Image Transport System\\) format',
                 '^ *COMMENT   Astronomy and Astrophysics, volume 376')

        for (string in reserved) {
            idx <- grep(string, header, ignore.case=TRUE)
            if (length(idx) > 0)  header <- header[-idx]
        }
        cimages <- c(cimages, header)
    }

    # Return complete and (if needed, merged) header
    closeHdr(cimages)
}
