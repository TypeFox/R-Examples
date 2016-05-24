readFITSheader <-
function (zz, maxLines = 5000, fixHdr = 'none')
{
### Function gets FITS header
###
### Takes:
  ## File handle: zz
  ## Maximum number of header lines to read: maxLines
### Returns:
  ## Character vector with header
###
### Refs: http://fits.gsfc.nasa.gov/
###       Hanisch et al., Astr. Ap. 376, 359-380 (2001)
###
### Modified version of readFITSheader.R from 3/21/08 through 9/25/2010
### A. Harris, Univ. MD Astronomy, 12/30/12
###

    ## Each header is a set of 36 80-column card images
    num <- 36      # number of card images in block
    cols <- 80     # columns per card
    maxHdrs <- maxLines/num  # max number of header units to read
    header <- 'dummy start line'
    image <- character(num)       # vector for card images
    start <- seq(1, 2880, by=80)  # start character index for images

    for (i in 1:maxHdrs) {
        # Read 36*80 = 2880 characters at a time from file header, with
        # different modes to accomodate standard header problems
        switch(pmatch(fixHdr[1], c('none', 'remove', 'substitute'), nomatch=4),
               { # header is ok
                    inpString <- readChar(zz, 2880)
                    if (nchar(inpString) != 2880) {
                        txt <- paste('*** Header problem:', nchar(inpString),
                        'characters instead of 2880; try option fixHdr *** \n')
                        close(zz)
                        stop(txt)
                    }
                }, { # substitute blanks for non-printing characters from header
                   nbytes <- 2880
                   inpBin <- raw()
                   while (nbytes != 0) {
                       inpBin <- c(inpBin, readBin(zz, 'raw', nbytes))
                       idxBad <- which(inpBin <= 0x1F | inpBin == 0x7F)
                       nbytes <- length(idxBad)
                       if (nbytes > 0) {
                           inpBin <- inpBin[-idxBad]
                           txt <- paste('*** Removed', length(idxBad),
                                    'non-printing characters in header ***\n')
                           cat(txt)
                       }
                       inpString <- rawToChar(inpBin)
                   }
                }, { # substitute spaces for non-printing characters
                    inpBin <- readBin(zz, 'raw', 2880)
                    idxBad <- which(inpBin <= 0x1F | inpBin == 0x7F)
                    nbytes <- length(idxBad)
                    if (nbytes > 0) {
                        inpBin[idxBad] <- as.raw(0x20)  # substitute space
                        txt <- paste('*** Substituted', nbytes,
                                     'spaces for non-printing characters ***\n')
                        cat(txt)
                    }
                    inpString <- rawToChar(inpBin)
                }, {
                    close(zz)
                    stop('*** Invalid fixHdr option ***\n')
                }
        )

        # Break header data into card images,
        for (j in 1:num) {
            image[j] <- substr(inpString, start[j], start[j]+79)
        }

        # Look for END card; if found, clean up and return header,
        # else add images to header
        idx <- grep('^ *END +', image, ignore.case=TRUE)
        if (length(idx) > 0) {
            image <- image[-(idx:num)]  # trim END and trailing blank cards
            header <- c(header, image)  # update header
            header <- header[-1]        # remove initial dummy line
            return(header)              # return
        } else {
            header <- c(header, image)  # update header
        }
    }

    # Return on error -- no END
    stop("Haven't found END in header after ", maxLines,
         " header lines")
}

