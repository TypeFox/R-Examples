read.fitskey = function(key, file, hdu = 1, comments = FALSE, strip = c(" ","'"," "), maxlines = 50000){
    
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
    
    # locate keys
    out = {}
    keys = hdr[,"key"]
    for(i in 1:length(key)){
        if(key[i]%in%keys){
            k = as.vector(hdr[which(keys==key[i]),"value"])
        }else{
            k = NA
        }
        out = c(out,k)
    }
    
    # return keys
    return(out)
    
}

