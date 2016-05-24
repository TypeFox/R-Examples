put.fitskey = function(key, value, hdr){
    
    # create key/value vectors
    k = rep(NA, max(length(key),length(value)))
    v = rep(NA, max(length(key),length(value)))
    k[1:length(key)] = key
    v[1:length(value)] = value
    
    # update pre-existing keys
    if(any(k %in% hdr[,"key"])){
        
        # determine order of vectors
        inhdr = which(hdr[,"key"] %in% k)
        inkey = which(k %in% hdr[,"key"])
        matched = inhdr[match(hdr[inhdr,"key"],k[inkey])]
        
        # update header
        hdr[matched,"value"] = v[inkey]
        
    }
    
    # create new keys
    if(any(!k %in% hdr[,"key"])){
        
        # determine new keys
        outhdr = which(!k %in% hdr[,"key"])
        
        # append to header
        hdr = rbind(hdr, c(key=k[outhdr], value=v[outhdr], comment=""))
        
    }
    
    # return results
    return(hdr)
    
}

