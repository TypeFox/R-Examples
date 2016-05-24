closeHdr <-
function(headerName) {
# Function adds end statement and closes FITS header
#
# Returns string array vector to write to FITS file
#
# A. Harris 2012.10.13
    headerName <- c(headerName, sprintf('END%-77s', ' '))
    tmp <- paste(headerName, collapse='')
    len <- nchar(tmp)
    nch <- 36*80 # card images times characters
    nblocks <- ceiling(len/nch)
    nfill <- nblocks*nch - len
    tmp <- paste(c(tmp, sprintf('%*s', nfill, ' ')), collapse='')

    out <- character(nblocks)
    sta <- 1
    sto <- nch
    for (i in 1:nblocks) {
        out[i] <- substr(tmp, sta, sto)
        sta <- sta + nch
        sto <- sto + nch
    }
    out
}
