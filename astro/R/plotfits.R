plotfits = function(input, hdu = 1, func = "atan", slide = c(0,0,0), scale = c(500,300,100), locut = 0, hicut = pi/2, invert = FALSE, method = 1, type = "x11", width = 5, height = 5, units = "in", res = 300, cen = c(NA,NA), xdim = NA, ydim = NA, file = "image.png"){
    
    #hdu = 1; func = "atan"; slide = c(0,0,0), scale = c(500,300,100); locut = 0; hicut = pi/2; invert = FALSE; method = 1; type = "x11"; width = 5; height = 5; units = "in"; res = 300; cen = c(NA,NA); xdim = NA; ydim = NA; file = "image.png"
    
    # format
    if(typeof(input)=="double"){
        input = list(input)
    }
    
    # bitmap types
    btypes = c("bit", "bitcmyk", "bitrgb", "jpeg", "jpeggray", "pbm", "pbmraw", "pcx16", "pcx24b", "pcx256", "pcxcmyk", "pcxgray", "pcxmono", "pdfwrite", "pgm", "pgmraw", "pgnm", "pgnmraw", "pkm", "pkmraw", "png16", "png16m", "png256", "png48", "pngalpha", "pnggray", "pngmono", "pnm", "pnmraw", "ppm", "ppmraw", "psgray", "psmono", "psrgb", "tiff12nc", "tiff24nc", "tiffcrle", "tiffg3", "tiffg32d", "tiffg4", "tifflzw", "tiffpack")
    
    # check inputs match
    if(length(slide)!=length(input)){slide=rep(slide[1],length(input))}
    if(length(scale)!=length(input)){scale=rep(scale[1],length(input))}
    if(length(locut)!=length(input)){locut=rep(locut[1],length(input))}
    if(length(hicut)!=length(input)){hicut=rep(hicut[1],length(input))}
    
    # image process function
    imageprocess = function(imdat, func, slide, scale, locut, hicut, cen, xdim, ydim){
        
        # read in data and scale
        if(typeof(imdat)=="character"){
            imdat = read.fitsim(imdat, hdu=hdu)
        }
        # centre
        if(is.na(cen[1])){
            cen = (dim(imdat)+1)/2
        }
        # xdim
        if(is.na(xdim)){
            xdim = dim(imdat)[1]
        }
        # ydim
        if(is.na(ydim)){
            ydim = dim(imdat)[2]
        }
        # populate fake data matrix
        fake = matrix(0,xdim,ydim)
        fcen = (dim(fake)+1)/2
        # trim imdat if necessary
        if(dim(imdat)[1]>dim(fake)[1]){
            imdat = imdat[(cen[1]-((xdim-1)/2)):(cen[1]+((xdim-1)/2)),]
            cen[1] = ((xdim+1)/2)
        }
        if(dim(imdat)[2]>dim(fake)[2]){
            imdat = imdat[,(cen[2]-((ydim-1)/2)):(cen[2]+((ydim-1)/2))]
            cen[2] = ((ydim+1)/2)
        }
        foff = fcen-cen
        xoff = (1:(dim(imdat)[1]))+(foff[1])
        yoff = (1:(dim(imdat)[2]))+(foff[2])
        fake[xoff,yoff] = imdat
        imdat = fake
        # apply slide
        imdat = imdat-slide
        # apply function
        if(func=="atan"){
            imdat = imdat*(1/scale)
            imdat = atan(imdat)
        }else if(func=="log"){
            imdat = imdat+500
            imdat = imdat*(1/scale)
            imdat[imdat<=0.001] = 0.001
            imdat = log10(imdat)
        }else{
            imdat = imdat*(1/scale)
        }
        # apply cut limits and scale from 0-1
        imdat[imdat<locut] = locut
        imdat[imdat>hicut] = hicut
        imdat = imdat-locut
        imdat = imdat/(hicut-locut)
        # return image data
        return(imdat)
        
    }
    
    # image manipulation and master matrix
    imblock = {}
    for(i in 1:length(input)){
        # generate imdat
        out = imageprocess(imdat=input[[i]], func=func, slide=slide[i], scale=scale[i], locut=locut[i], hicut=hicut[i], cen=cen, xdim=xdim, ydim=ydim)
        imblock = c(imblock, list(out))
    }
    mat=matrix(1:length(imblock[[1]]),dim(imblock[[1]])[1],dim(imblock[[1]])[2])
    
    # colourify
    if(method==1){
        # RGB (simple addition): scale imdat from 1:256
        if(length(imblock)==1){
            imblock = c(imblock,list(imblock[[1]],imblock[[1]]))
        }
        # scale 0:1 -> 1:256
        red = (imblock[[1]]*255)+1
        green = (imblock[[2]]*255)+1
        blue = (imblock[[3]]*255)+1
        # colours
        clim = c(0,1)
        if(invert){
            clim = c(1,0)
        }
        # return colour
        col=rgb(seq(clim[1],clim[2],len=256)[red],seq(clim[1],clim[2],len=256)[green],seq(clim[1],clim[2],len=256)[blue])
    }else{
        # HSV (complex colour wheel): variable saturation - under construction
        if(length(imblock)==1){
            imblock = c(imblock,list(0,0))
        }
        # polar -> cartesian
        xred = imblock[[1]]*cos((0)*(2*pi))
        yred = imblock[[1]]*sin((0)*(2*pi))
        xgreen = imblock[[2]]*cos((1/3)*(2*pi))
        ygreen = imblock[[2]]*sin((1/3)*(2*pi))
        xblue = imblock[[3]]*cos((2/3)*(2*pi))
        yblue = imblock[[3]]*sin((2/3)*(2*pi))
        # find mean positions
        xnew = (xred+xgreen+xblue)/3
        ynew = (yred+ygreen+yblue)/3
        # cartesian -> polar
        tweak = 0
        if(invert){tweak=180}
        hue = (((atan2(ynew,xnew)*(180/pi))+tweak)%%360)/360
        sat = sqrt((xnew^2) + (ynew^2)) * 3
        # return colour
        col = hsv(h=hue,s=sat)
    }
    
    # plotting dimensions
    if(is.na(width) | is.na(height) | toupper(width)=="NA" | toupper(height)=="NA"){
        width = dim(imblock[[1]])[2]
        height = dim(imblock[[1]])[1]
        units = "px"
    }
    
    # plot
    if(type=="png"){
        png(width=width, height=height, units=units, res=res, filename=file)
    }else if(type%in%btypes){
        bitmap(type=type, width=width, height=height, units=units, res=res, file=file)
    }else if(type=="eps"){
        setEPS()
        postscript(width=width, height=height, file=file)
    }else if(type=="pdf"){
        pdf(width=width, height=height, file=file)
    }
    if(type!="dat"){
        par(mar=c(0,0,0,0), mgp=c(3,1,0), col="white", bg="black")
        image(mat, col=col, axes=FALSE)
    }
    if(type!="x11" & type!="X11" & type!="dat"){
        graphics.off()
    }
    # return data frame and colours
    if(type=="dat"){
        return(list(mat=mat,col=col))
    }
    
}

