write.fitshdr = function(hdr, file, hdu = 1){
    
    # read in original file
    x = read.fits(file, hdu=0, comments=TRUE, strip=c(" ","'"," "))
    
    # comments?
    if(length(hdr[1,])==2){
        hdr = cbind(hdr,comment="")
    }
    
    # update x
    x$hdr[[hdu]] = hdr
    
    # write new file
    write.fits(x, file = file, type = "auto")
    
}

