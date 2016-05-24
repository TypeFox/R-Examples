read.fits = function(file, hdu = 0, comments = TRUE, strip = c(" ","'"," "), maxlines = 50000, xlo = NA, xhi = NA, ylo = NA, yhi = NA){
    
    # open file
    fcon = file(file, "rb")
    
    # read primary header
    hdr = .read.fits.hdr(fcon, comments=comments, strip=strip, maxlines=maxlines)
    
    # read primary data unit
    du = .read.fits.image(fcon, hdr, xlo, xhi, ylo, yhi)
    
    # correct hdr for WCS information if present and if cutting image
    wcspix = as.numeric(get.fitskey(key=c("CRPIX1","CRPIX2"), hdr=hdr))
    if(!is.na(xlo)){
        if(!is.na(wcspix[1])){
            wcspix[1] = wcspix[1] - xlo + 1
            hdr = put.fitskey(key="CRPIX1", value=wcspix[1], hdr=hdr)
        }
    }
    if(!is.na(ylo)){
        if(!is.na(wcspix[2])){
            wcspix[2] = wcspix[2] - ylo + 1
            hdr = put.fitskey(key="CRPIX2", value=wcspix[2], hdr=hdr)
        }
    }
    
    # master files
    hlist = list(hdr)
    dlist = list(du)
    
    # secondaries
    xhdu = TRUE
    nhdu = 1
    if(hdu==0){thdu=Inf}else{thdu=hdu}
    while(xhdu & thdu>nhdu){
        
        # read header
        hdr = .read.fits.hdr(fcon, comments=comments, strip=strip, maxlines=maxlines)
        
        # check hdu exists
        if(length(hdr)==1){
            if(!hdr[1]){
                xhdu = FALSE
            }
        }
        
        # read data unit
        if(xhdu){
            
            # read fits table/image
            if(length(grep("BINTABLE",hdr[which(hdr[,"key"]=="XTENSION"),"value"]))>0){
                du = .read.fits.table(fcon, hdr)
            }else{
                du = .read.fits.image(fcon, hdr, xlo, xhi, ylo, yhi)
            }
            
            # correct hdr for WCS information if present and if cutting image
            wcspix = as.numeric(get.fitskey(key=c("CRPIX1","CRPIX2"), hdr=hdr))
            if(!is.na(xlo)){
                if(!is.na(wcspix[1])){
                    wcspix[1] = wcspix[1] - xlo + 1
                    hdr = put.fitskey(key="CRPIX1", value=wcspix[1], hdr=hdr)
                }
            }
            if(!is.na(ylo)){
                if(!is.na(wcspix[2])){
                    wcspix[2] = wcspix[2] - ylo + 1
                    hdr = put.fitskey(key="CRPIX2", value=wcspix[2], hdr=hdr)
                }
            }
            
            # add to masters
            hlist = c(hlist, list(hdr))
            dlist = c(dlist, list(du))
            
            # header count
            nhdu = nhdu+1
            
        }
        
    }
    
    # close file
    close(fcon)
    
    # choose hdu
    if(hdu>0){
        hlist = list(hlist[[hdu]])
        dlist = list(dlist[[hdu]])
    }
    
    # return
    return(list(hdr=hlist,dat=dlist))
    
}

.read.fits.hdr = function(fcon, comments=TRUE, strip=c(" ","'"," "), maxlines=50000){
    
    # setup
    hdr = character()
    i = 0
    maxheadroom = maxlines/36
    foundend = FALSE
    
    ## FITS headers are always multiples of 36 rows (*80 columns = 2880 character chunks)
    while(!foundend){
        rawhead = readChar(fcon, 36*80)
        if(length(rawhead)>0){
            #rawhead = readBin(fcon, what="character")
            temp = .parse.fits.hdr(rawhead, comments=comments, strip=strip)
            hdr = rbind(hdr, temp)
            if(any(hdr[,"key"]=="END")){
                end = which(hdr[,"key"]=="END")
                hdr = hdr[-(end:length(hdr[,"key"])),]
                foundend = TRUE
            }
            if(i > maxheadroom){
                stop("Haven't found END in header after ", maxlines, " header lines")
            }
            i = i + 1
        }else{
            foundend = TRUE
            hdr = FALSE
        }
    }
    
    # return header
    return(hdr)
    
}

.parse.fits.hdr = function(rawhead, comments=TRUE, strip=c(" ","'"," ")){
    
    # create and parse header rows
    key = {}
    value = {}
    comment = {}
    for(i in 1:(nchar(rawhead)/80)){
    
        # row
        row = substr(rawhead,((80*(i-1))+1),(80*i))
        
        # key
        k = paste(strsplit(substr(row,1,8)," +")[[1]],collapse="",sep="")
        if(k=="HIERARCH"){
            charend = (min(which(strsplit(row,"")[[1]]=="="))-2)
            k = strip(paste(strsplit(row,"")[[1]][10:charend],collapse="",sep=""), strip=strip)
            leftover = substr(row,charend+3,80)
        }else{
            leftover = substr(row,10,80)
        }
        key = c(key, k)
        
        # value(?)
        if(k=="COMMENT" | k=="HISTORY"){
            v = ""
        }else{
            poss = strsplit(leftover,"/")[[1]]
            while(((length(grep("'",strsplit(poss[1],"")[[1]]))%%2) != 0) & (length(poss)>2)){
                poss[1] = paste(poss[1],"/",poss[2],collapse="",sep="")
                poss = poss[-2]
            }
            temp = strip(poss[1], strip=strip)
            v = paste(temp,collapse=" ")
        }
        value = c(value,v)
        
        # comment(?)
        if(comments){
            if(k=="COMMENT" | k=="HISTORY"){
                m = strip(substr(row,9,80), strip=" ")
            }else{
                m = strip(paste(poss[2:length(poss)],collapse="/"), strip=" ")
            }
            comment = c(comment,m)
        }
    }
    
    # collect results
    if(comments){
        hdr = cbind(key,value,comment)
    }else{
        hdr = cbind(key,value)
    }
    return(hdr)
    
}

.read.fits.image = function(fcon, hdr, xlo, xhi, ylo, yhi){
    
    # choose data precision type
    switch(hdr[which(hdr[,"key"]=="BITPIX"),"value"], "-64" = {
        bsize = 8           # 64-bit float (double)
        btype = numeric()
        bsign = TRUE
    }, "-32" = {
        bsize = 4           # 32-bit float (single)
        btype = numeric()
        bsign = TRUE
    }, "32" = {
        bsize = 4           # 32-bit signed int (double)
        btype = integer()
        bsign = TRUE
    }, "16" = {
        bsize = 2           # 16-bit signed int (single)
        btype = integer()
        bsign = TRUE
    }, "8" = {
        bsize = 1           # 8-bit unsigned int
        btype = integer()
        bsign = FALSE
    }, stop("Unknown BITPIX request in .read.fits.image"))
    
    # put array information into vectors
    naxis = as.numeric(hdr[which(hdr[,"key"]=="NAXIS"),"value"])
    
    # if image information contained
    if(naxis!=0){
        naxisn = integer(naxis)
        records = 1
        for(i in 1:naxis){
            temp = as.numeric(hdr[which(hdr[,"key"]==paste("NAXIS",i,sep="")),"value"])
            records = records * temp
            naxisn[i] = temp
        }
        
        # read image data
        subregion = FALSE
        if(!is.na(xlo) & !is.na(xhi) & !is.na(ylo) & !is.na(yhi)){
            xcols = xlo:xhi
            ycols = ylo:yhi
            if(prod(naxisn) != prod(c(length(xcols),length(ycols)))){
                subregion = TRUE
            }
        }
        if(subregion){
            start = (naxisn[1] * (ylo - 1)) + (xlo - 1)
            if(start != 0){seek(fcon, where = start*bsize, origin = "current")}
            dat = {}
            for(row in 1:length(ycols)){
                temp = readBin(fcon, what = btype, n = length(xcols), size = bsize, signed = bsign, endian = "big")
                dat = cbind(dat, temp)
                if(row != yhi){
                    seek(fcon, where = (naxisn[1] - length(xcols)) * bsize, origin = "current")
                }
            }
            dimnames(dat)[[2]] = NULL
            finish = (naxisn[1] * (naxisn[2] - yhi)) + (naxisn[1] - xhi)
            if(finish != 0){seek(fcon, where = finish*bsize, origin = "current")}
        }else{
            dat = array(readBin(fcon, what = btype, n = records, size = bsize, signed = bsign, endian = "big"), dim = naxisn)
        }
        
        # finish reading block to allow reading next hdu (if exists)
        nbyte = records * bsize
        nbyte = ifelse(nbyte%%2880 == 0, 0, 2880 - nbyte%%2880)
        temp = readBin(fcon, what = 'raw', n = nbyte)
        
        # scale image to physical units if needed
        temp = hdr[which(hdr[,"key"] == "BUNIT"),"value"]
        BUNIT = ifelse(length(temp) != 1, "", temp)
        temp = hdr[which(hdr[,"key"] == "BSCALE"),"value"]
        BSCALE = ifelse(length(temp) != 1, 1, as.numeric(temp))
        temp = hdr[which(hdr[,"key"] == "BZERO"),"value"]
        BZERO = ifelse(length(temp) != 1, 0, as.numeric(temp))
        if(BSCALE != 1 || BZERO != 0){
            dat = dat * BSCALE + BZERO
        }
        
    }else{
        dat = NULL
    }
        
    # return image data
    return(dat)
    
}

.read.fits.table = function(fcon, hdr){
    
    # setup
    naxis1 = as.numeric(hdr[which(hdr[,"key"]=="NAXIS1"),"value"])      # bytes per row
    naxis2 = as.numeric(hdr[which(hdr[,"key"]=="NAXIS2"),"value"])      # rows
    tfields = as.numeric(hdr[which(hdr[,"key"]=="TFIELDS"),"value"])    # columns
    
    # Get additional memory allocation information
    pcount = as.numeric(hdr[which(hdr[,"key"]=="PCOUNT"),"value"])      # size of special data area
    gcount = as.numeric(hdr[which(hdr[,"key"]=="GCOUNT"),"value"])      # one data group
    if(length(pcount)!=1){pcount = 0}
    if(length(gcount)!=1){gcount = 1}
    if(pcount != 0){warning("pcount must be 0 in a bintable")}
    if(gcount != 1){warning("gcount must be 1 in a bintable")}
    
    # Put array information into vectors
    # used: TTYPE, TFORM, TUNIT, TCOMM, TSCAL, TZERO
    # unused: TUCD, TNULL, TDISP, THEAP, TDIM
    ttypen = tformn = tunitn = tcommn = character(tfields)
    tscaln = tzeron = numeric(tfields)
    for(i in 1:tfields){
        # TTYPE - column name
        ttype = hdr[which(hdr[,"key"]==paste("TTYPE",i,sep="")),"value"]
        if(length(ttype)>0){
            ttype = strip(strip(ttype,strip="'")," ")
        }else{
            ttype = ""
        }
        ttypen[i] = ttype
        # TFORM - format type
        tform = hdr[which(hdr[,"key"]==paste("TFORM",i,sep="")),"value"]
        if(length(tform)>0){
            tform = strip(strip(tform,strip="'")," ")
        }else{
            tform = ""
        }
        tformn[i] = tform
        # TUNIT - units
        tunit = hdr[which(hdr[,"key"]==paste("TUNIT",i,sep="")),"value"]
        if(length(tunit)>0){
            tunit = strip(strip(tunit,strip="'")," ")
        }else{
            tunit = ""
        }
        tunitn[i] = tunit
        # TCOMM - column descriptions
        tcomm = hdr[which(hdr[,"key"]==paste("TCOMM",i,sep="")),"value"]
        if(length(tcomm)>0){
            tcomm = strip(strip(tcomm,strip="'")," ")
        }else{
            tcomm = ""
        }
        tcommn[i] = tcomm
        # TSCAL
        tscal = hdr[which(hdr[,"key"]==paste("TSCAL",i,sep="")),"value"]
        if(length(tscal)>0){
            tscal = strip(strip(tscal,strip="'")," ")
        }else{
            tscal = 1
        }
        tscaln[i] = tscal
        # TZERO
        tzero = hdr[which(hdr[,"key"]==paste("TZERO",i,sep="")),"value"]
        if(length(tzero)>0){
            tzero = strip(strip(tzero,strip="'")," ")
        }else{
            tzero = 0
        }
        tzeron[i] = tzero
    }
    
    # Work out formats and set up storage for each column.
    bsize = integer(tfields)
    btype = integer(tfields)
    bsign = logical(tfields)
    mult = integer(tfields)
    col = vector("list", tfields)
    tformcharn = {}
    for(i in 1:tfields){
        
        # TFORM handler
        tform = strsplit(tformn[i],"")[[1]]
        if(any(tform=="(")){
            brac = which(tform=="(")
            tform = tform[1:(brac-1)]
        }
        if(any(tform=="P")){
            pnum = which(tform=="P")
            tform = tform[-pnum]
        }
        tformchar = tform[length(tform)]
        tform = tform[-which(tform%in%tformchar)]
        tformcharn = c(tformcharn, tformchar)
        tformnumr = as.numeric(paste(tform, collapse=""))
        if(is.na(tformnumr)){tformnumr = 1}
        if(length(tformnumr)==0){tformnumr = 1}
        mult[i] = tformnumr
        
        # format; btype: 1 character, 2 logical, 3 integer, 4 numeric, 5 complex
        form = tolower(tformchar)
        switch(form, l = {
            bsize[i] = 1    # 1-byte/8-bit logical/boolean 
            btype[i] = 2
            bsign[i] = FALSE
        }, b = {
            bsize[i] = 1    # 1-byte/8-bit unsigned int
            btype[i] = 3
            bsign[i] = FALSE
        }, i = {
            bsize[i] = 2    # 2-byte/16-bit signed int (short)
            btype[i] = 3
            bsign[i] = TRUE
        }, j = {
            bsize[i] = 4    # 4-byte/32-bit signed int (long)
            btype[i] = 3
            bsign[i] = TRUE
        }, k = {
            bsize[i] = 8    # 8-byte/64-bit signed int (long long)
            btype[i] = 3    # does not entirely work due to an R 64-bit issue?
            bsign[i] = TRUE
        }, a = {
            bsize[i] = 1    # 1-byte/8-bit character
            btype[i] = 1
            bsign[i] = FALSE
        }, e = {
            bsize[i] = 4    # 4-byte/32-bit float (single)
            btype[i] = 4
            bsign[i] = TRUE
        }, d = {
            bsize[i] = 8    # 8-byte/64-bit float (double)
            btype[i] = 4
            bsign[i] = TRUE
        }, stop("X, C, M, P FITS format codes not yet implemented \n"))
        
        ## Set up storage arrays: rows = number of table rows, columns =
        ## multiplier (depth) of the cells in each column.  Characters are an
        ## exception since they return strings.
        if(btype[i] == 1){
            col[[i]] = array("", dim = c(naxis2, 1))
        }else{
            col[[i]] = array(switch(btype[i], "", FALSE, NA, NA, NA), dim = c(naxis2, mult[i]))
        }
        
    }
    
    # Read data, row by row
    for(i in 1:naxis2){
        for(j in 1:tfields){
            if(btype[j] <= 2){  # character reads
                col[[j]][i,] = readChar(fcon, nchars = mult[j])
            }else{
                what = switch(btype[j], character(), logical(), integer(), numeric(), complex())
                col[[j]][i,] = readBin(fcon, what = what, n = mult[j], size = bsize[j], signed = bsign[j], endian = "big")
            }
        }
    }
    
#    # EXPERIMENTAL SECTION TO DEAL WITH COMPRESSED FITS IMAGES
#    # heap area
#    if(pcount>0){
#        temp = readBin(fcon, 'raw', pcount)
#    }
    
    # finish reading block to allow reading next hdu (if exists)
    nbyte = (naxis1 * naxis2)
    nbyte = ifelse(nbyte%%2880 == 0, 0, 2880 - nbyte%%2880)
    temp = readBin(fcon, 'numeric', nbyte)
    
    ## Clean up before returning
    for(i in 1:tfields){
        ## Apply scaling and offset where appropriate
        if(btype[i] >= 3 && (tscaln[i] != 1 || tzeron[i] != 0)){
            col[[i]] = (col[[i]]) * tscaln[i] + tzeron[i]
        }
        ## Convert 1D arrays to vectors for easier plotting
        if(nrow(col[[i]]) == 1 || ncol(col[[i]]) == 1 || btype[i] <= 2){
            col[[i]] = as.vector(col[[i]])
        }
#        ## Terminate character strings that end in ascii 0
# commented out 2 Mar 2012, causing problems - splits character strings up into multiple strings if it contains a zero.
#        if(btype[i] <= 2){
#            lcol = length(col[[i]])
#            txttmp = character(lcol)
#            for(j in 1:lcol){
#                txttmp[j] = strsplit((col[[i]][j]), 0)
#            }
#            col[[i]] = unlist(txttmp)
#        }
    }
    
    # format outputs
    table = {}
    for(i in 1:tfields){
        table = cbind(table,unlist(col[[i]]))
    }
    colnames(table)=ttypen
    meta = cbind(name=ttypen,units=tunitn,description=tcommn)
    
    # return outputs
    return(list(meta=meta,table=table))
    
}

