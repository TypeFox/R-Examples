write.fits = function(x, file = "star.fits", type = "single", hdu = 0){
    
    # split x into hdr/dat
    if(!is.list(x)){x=list(x)}
    if(!is.null(x$hdr)){
        hdr = x$hdr
        dat = x$dat
    }else{
        hdr = rep(NA,length(x))
        dat = x
    }
    if(!is.list(hdr)){hdr = as.list(hdr)}
    if(hdu > 0){
        hdr = list(hdr[[hdu]])
        dat = list(dat[[hdu]])
    }
    
    # open file
    fcon = file(file, "wb")
    
    # loop over each hdu in dat
    for(i in 1:length(dat)){
        
        # sub-setup
        if(is.null(dat[[i]])){
            d = dat[[i]]
        }else{
            d = as.array(dat[[i]])
        }
        h = hdr[[i]]
        naxis = length(dim(d))
        
        type = tolower(substr(type, 1, 1))
        
        if(!is.na(h[1]) & type=="a"){
            # choose data precision type from header
            bitpix = as.numeric(h[which(h[,"key"]=="BITPIX"),"value"])
            switch(as.character(bitpix), "-64" = {
                size = 8    # 64-bit float (double)
            }, "-32" = {
                size = 4    # 32-bit float (single)
            }, "32" = {
                size = 4    # 32-bit signed int (double)
            }, "16" = {
                size = 2    # 16-bit signed int (single)
            }, "8" = {
                size = 1    # 8-bit unsigned int
            }, stop("Unknown BITPIX request in write.fits"))
        }else{
            # Choose integer/float; byte/single/double precision based on data
            if(is.null(d)){
                bitpix = 8
                size = 0
            }else if(is.integer(d)){
                switch(type, b = {
                    bitpix = 8
                    size = 1
                }, s = {
                    bitpix = 16
                    size = 2
                }, d = {
                    bitpix = 32
                    size = 4
                }, stop("Unrecognized data type"))
            }else{
                switch(type, s = {
                    bitpix = -32
                    size = 4
                }, d = {
                    bitpix = -64
                    size = 8
                }, stop("Unrecognized data type"))
            }
        }
        
        # dummy header if none provided/check header provided
        if(i==1){phdu = TRUE}else{phdu = FALSE}
        if(is.na(h[1])){
            h = .dummy.fits.hdr(phdu=phdu, bitpix=bitpix, naxis=naxis, naxisn=dim(d))
        }else{
            h = .check.fits.hdr(hdr=h, phdu=phdu, bitpix=bitpix, naxis=naxis, naxisn=dim(d))
        }
        
        # make fits header
        h = .make.fits.hdr(h)
        
        # write primary header
        writeChar(h, fcon, eos = NULL)
        
        if(!is.null(d)){
            
            # write primary data unit
            writeBin64(as.vector(d), fcon, size = size, endian = "big")
            
            # pad rest of record with zeros
            pad = raw(2880 - (length(as.vector(d)) * size)%%2880)
            writeBin(pad, fcon, endian = "big")
            
        }
        
    }
    
    # close file
    close(fcon)
    
}

#write.fitsim = function(){
#    
#}

#write.fitstab = function(){
#    
#}

.dummy.fits.hdr = function(phdu=TRUE, bitpix, naxis, naxisn){
    
    # naxis
    naxisnames = {}
    naxiscomm = {}
    for(i in 1:naxis){
        naxisnames = c(naxisnames, paste("NAXIS",i,sep=""))
        naxiscomm = c(naxiscomm, paste("length of data axis",i))
    }
    
    # setup base
    key = c("BITPIX", "NAXIS", naxisnames)
    value = c(bitpix, naxis, naxisn)
    comment = c("number of bits per data pixel", "number of data axes", naxiscomm)
    
    # primary/secondary header
    if(phdu){
        
        key = c("SIMPLE", key, "EXTEND", "COMMENT", "COMMENT")
        value = c("T", value, "T", "", "")
        comment = c("file does conform to FITS standard", comment, "FITS dataset may contain extensions", "FITS (Flexible Image Transport System) format is defined in 'Astronomy", "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H")
        
    }else{
        
        key = c("XTENSION", key, "PCOUNT", "GCOUNT")
        value = c("IMAGE", value, "0", "1")
        comment = c("IMAGE extension", comment, "number of random group parameters", "number of random groups")
        
    }
    
    # create and return header
    hdr = cbind(key,value,comment)
    return(hdr)
    
}

.check.fits.hdr = function(hdr, phdu, bitpix, naxis, naxisn){
    
    # comments
    if(!"comment"%in%colnames(hdr)){
        c = rep("",length(hdr[,1]))
        hdr = cbind(hdr,comment=c)
    }
    
    # naxis
    naxisnames = {}
    naxiscomm = {}
    if(naxis>0){
        for(i in 1:naxis){
            naxisnames = c(naxisnames, paste("NAXIS",i,sep=""))
            naxiscomm = c(naxiscomm, paste("length of data axis",i))
        }
    }
    
    # primary/secondary
    if(phdu){
        extras = c("SIMPLE","EXTEND")
        xvals = c("T","T")
        xcoms = c("file does conform to FITS standard","FITS dataset may contain extensions")
        remove = c("XTENSION","PCOUNT","GCOUNT")
    }else{
        extras = c("XTENSION","PCOUNT","GCOUNT")
        xvals = c("IMAGE","0","1")
        xcoms = c("IMAGE extension","number of random group parameters","number of random groups")
        remove = c("SIMPLE","EXTEND")
    }
    
    # check columns and values
    colcheck = c(extras, "BITPIX", "NAXIS", naxisnames)
    colvals = c(xvals, bitpix, naxis, naxisn)
    colcoms = c(xcoms, "number of bits per data pixel", "number of data axes", naxiscomm)
    
    for(i in length(colcheck):1){
        if(colcheck[i]%in%hdr[,"key"]){
            col = which(hdr[,"key"]==colcheck[i])
            if(as.character(hdr[col,"value"])!=as.character(colvals[i])){
                hdr[col,"value"] = colvals[i]
            }
        }else{
            hdr = rbind(c(colcheck[i],colvals[i],colcoms[i]),hdr)
        }
    }
    for(j in 1:length(remove)){
        if(remove[j]%in%hdr[,"key"]){
            col = which(hdr[,"key"]==remove[j])
            hdr = hdr[-col,]
        }
    }
    
    # reorder main keys
    if(phdu){
        corder = c("SIMPLE","BITPIX","NAXIS",naxisnames,"EXTEND")
    }else{
        corder = c("XTENSION","BITPIX","NAXIS",naxisnames,"PCOUNT","GCOUNT")
    }
    temp = {}
    hdr = rbind(hdr,c(0,0,0),c(0,0,0))
    for(k in 1:length(corder)){
        if(corder[k]%in%hdr[,"key"]){
            col = which(hdr[,"key"]==corder[k])
            temp = rbind(temp,hdr[col,])
            hdr = hdr[-col,]
        }
    }
    hdr = rbind(temp,hdr)
    bad = (length(hdr[,1])-1):(length(hdr[,1]))
    hdr = hdr[-bad,]
    
    # return
    return(hdr)
    
}

.make.fits.hdr = function(hdr){
    
    # columns
    key = hdr[,"key"]
    value = hdr[,"value"]
    comment = hdr[,"comment"]
    
    # loop over each key
    rows = {}
    for(i in 1:length(key)){
    
        # special key?
        special = FALSE
        if(key[i]=="COMMENT" | key[i]=="HISTORY"){
            special = TRUE
        }
        
        # key
        key[i] = sprintf("%-8s", key[i])
        hier = FALSE
        if(nchar(key[i])>8){
            key[i] = paste("HIERARCH",key[i],"")
            hier = TRUE
        }
        
        # value
        if(special){
            value[i] = ""
        }else{
            # treat as a character?
            char = FALSE
            if(is.na(value[i])){
                char = TRUE
                value[i] = paste("'",sprintf("%-8s", value[i]),"'",sep="")
            }else if(is.na(suppressWarnings(as.numeric(value[i]))) & value[i]!="T" & value[i]!="F"){
                char = TRUE
                value[i] = paste("'",sprintf("%-8s", value[i]),"'",sep="")
            }
            # hierarch or normal?
            if(hier){
                pad = 33-nchar(key[i])-2-nchar(value[i])-3
                if(pad<0){pad=0}
                if(char){
                    value[i] = paste("= ", value[i],paste(rep(" ",pad),sep="",collapse=""), " / ", sep="", collapse="")
                }else{
                    value[i] = paste("= ", paste(rep(" ",pad),sep="",collapse=""),value[i], " / ", sep="", collapse="")
                }
            }else{
                if(char){
                    value[i] = paste("= ", sprintf("%-20s", value[i]), " / ", sep="")
                }else{
                    value[i] = paste("= ", sprintf("%20s", value[i]), " / ", sep="")
                }
            }
        }
        
        # update rows
        row = substr(sprintf("%-80s", paste(key[i], value[i], comment[i], sep="")),1,80)
        rows = c(rows,row)
        
    }
    
    # pad header to multiple of 2880 characters (36*80)
    hdr = paste(paste(rows,collapse="",sep=""),"END",sep="")
    pad = paste(rep(" ",(2880 - (nchar(hdr)%%2880))),collapse="",sep="")
    hdr = paste(hdr, pad, collapse="", sep="")
    
    return(hdr)
    
}

