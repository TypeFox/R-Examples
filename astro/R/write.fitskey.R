write.fitskey = function(key, value, file, comment = "", hdu = 1){
    
    # read in original file
    x = read.fits(file=file, hdu=0, comments=TRUE, strip=c(" ","'"," "))
    
    # check for key and add to hdr
    h = x$hdr[[hdu]]
    k = key
    v = m = character(length(k))
    if(length(value)>0){v[1:min(c(length(k),length(value)))] = value[1:min(c(length(k),length(value)))]}
    if(length(comment)>0){m[1:min(c(length(k),length(comment)))] = comment[1:min(c(length(k),length(comment)))]}
    for(i in 1:length(k)){
        if(k[i]%in%h[,"key"] & k[i]!="COMMENT" & k[i]!="HISTORY"){
            col = which(h[,"key"]==k[i])
            h[col,"value"] = v[i]
            h[col,"comment"] = m[i]
        }else{
            h = rbind(h,c(k[i],v[i],m[i]))
        }
    }
    x$hdr[[hdu]] = h
    
    # write out new file with updated header
    write.fits(x=x, file=file, type="auto")
    
}

