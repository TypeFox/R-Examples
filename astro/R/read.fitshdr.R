read.fitshdr = function(file, hdu = 1, comments = TRUE, strip = c(" ","'"," "), maxlines = 50000){
    
    if(hdu==1){
        
        # open file
        fcon = file(file, "rb")
        
        # read primary header
        hdr = .read.fits.hdr(fcon, comments=comments, strip=strip, maxlines=maxlines)
        
        # close file
        close(fcon)
        
    }else{
        
        # read and return header
        dat = read.fits(file, hdu=hdu, comments=comments, strip=strip, maxlines=maxlines)
        hdr = dat$hdr[[1]]
        
    }
    
    # return header
    return(hdr)
    
}

