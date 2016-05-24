read.fitstab = function(file, hdu=2, strip = c(" ","'"," "), maxlines = 50000){
    
    # read and return FITS table
    dat = read.fits(file, hdu=hdu, comments=FALSE, strip=strip, maxlines=maxlines)
    return(dat$dat[[1]]$table)
    
}

