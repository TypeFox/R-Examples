read.fitsim = function(file, hdu = 1, maxlines = 50000, xlo = NA, xhi = NA, ylo = NA, yhi = NA){
    
    # read and return FITS image
    dat = read.fits(file, hdu=hdu, comments=FALSE, maxlines=maxlines, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi)
    return(dat$dat[[1]])
    
}

