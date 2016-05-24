edges <- function(img, type="Laplacian", ltype=1, abs=FALSE){
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(img$compressed) img <- decompress.image(img)
  eimg <- switch(type,
                 "Laplacian" = laplacian(img$img,ltype),
                 "Sobel" = sobel(img$img),
                 "Robertcross" = robertcross(img$img), 
                 warning("Wrong type: must be Laplacian, Sobel or Robertcross"))
  deimg <- dim(eimg)
  if(abs)   eimg <- array(abs(eimg),deimg)
  if(length(deimg)==2){
   if((z <- diff(range(eimg))) > 0) eimg <- (eimg - min(eimg))/z
  } else {
     for( i in 1:deimg[3])
   if((z <- diff(range(eimg[,,i]))) > 0) eimg[,,i] <- (eimg[,,i] - min(eimg[,,i]))/z
   }
  invisible(make.image(eimg))
}

laplacian <- function(img,ltype=4) {
  dimg <- dim(img)
  ldim <- length(dimg)
  if (!(ldim %in% 2:3)) return (warning("Not an image"))
  lchannel <- dimg[1]*dimg[2]
  
  conv <- switch (as.character(ltype),
                  "1" = matrix(c(-1,-1,-1,-1,-1,
                                 -1,-1,-1,-1,-1,
                                 -1,-1,24,-1,-1,
                                 -1,-1,-1,-1,-1,
                                 -1,-1,-1,-1,-1), 5, 5),
                  "2" = matrix(c(0,-1,0,-1,4,-1,0,-1,0), 3, 3),
                  "3" = matrix(c(-1,-1,-1,-1,8,-1,-1,-1,-1), 3, 3),
                  matrix(c(1,-2,1,-2,4,-2,1,-2,1), 3, 3)
                  )
  
  eimg <- array(0,dim=dimg)
  if( ldim == 2 ) dim(eimg) <- dim(img) <- c(dimg, 1)
  for ( i in 1:dim(eimg)[3]) 
    eimg[,,i] <- .Fortran("convolve",
                          as.double(img[,,i]),
                          as.double(conv),
                          eimg=double(lchannel),
                          as.integer(dimg[2]),
                          as.integer(dimg[1]),
                          as.integer(dim(conv)[1]),
                          PACKAGE="adimpro")$eimg
  
  dim(eimg) <- dimg
  invisible(eimg)
}

sobel <- function(img) {
  dimg <- dim(img)
  ldim <- length(dimg)
  if (!(ldim %in% 2:3)) return (warning("Not an image"))
  lchannel <- dimg[1]*dimg[2]
  horizontal <- matrix(c(-1,0,1,
                         -2,0,2,
                         -1,0,1), 3, 3)
  vertical <- matrix(c(-1,-2,-1,
                       0,0,0,
                       1,2,1), 3, 3)
  
  ximg <- yimg <- array(0,dim=dimg)
  if( ldim == 2 ) dim(ximg) <- dim(yimg) <- dim(img) <- c(dimg, 1)
  for ( i in 1:dim(ximg)[3]) 
    ximg[,,i] <- .Fortran("convolve",
                          as.double(img[,,i]),
                          as.double(horizontal),
                          eimg=double(lchannel),
                          as.integer(dimg[2]),
                          as.integer(dimg[1]),
                          as.integer(dim(horizontal)[1]),
                          PACKAGE="adimpro")$eimg
  dim(ximg) <- dimg
  
  for ( i in 1:dim(yimg)[3]) 
    yimg[,,i] <- .Fortran("convolve",
                          as.double(img[,,i]),
                          as.double(vertical),
                          eimg=double(lchannel),
                          as.integer(dimg[2]),
                          as.integer(dimg[1]),
                          as.integer(dim(vertical)[1]),
                          PACKAGE="adimpro")$eimg
  dim(yimg) <- dimg
  
  invisible(sqrt(ximg^2 + yimg^2))
}

robertcross <- function(img) {
  dimg <- dim(img)
  ldim <- length(dimg)
  if (!(length(dimg) %in% 2:3)) return (warning("Not an image"))
  lchannel <- dimg[1]*dimg[2]
  horizontal <- matrix(c(1,0,
                         0,-1),
                       nrow=2,
                       ncol=2)
  vertical <- matrix(c(0,-1,
                       1,0),
                     nrow=2,
                     ncol=2)
  
  ximg <- yimg <- array(0,dim=dimg)
  if( ldim == 2 ) dim(ximg) <- dim(yimg) <- dim(img) <- c(dimg, 1)
  for ( i in 1:dim(ximg)[3]) 
    ximg[,,i] <- .Fortran("convolve",
                          as.double(img[,,i]),
                          as.double(horizontal),
                          eimg=double(lchannel),
                          as.integer(dimg[2]),
                          as.integer(dimg[1]),
                          as.integer(dim(horizontal)[1]),
                          PACKAGE="adimpro")$eimg
  dim(ximg) <- dimg
  
  for ( i in 1:dim(yimg)[3]) 
    yimg[,,i] <- .Fortran("convolve",
                          as.double(img[,,i]),
                          as.double(vertical),
                          eimg=double(lchannel),
                          as.integer(dimg[2]),
                          as.integer(dimg[1]),
                          as.integer(dim(vertical)[1]),
                          PACKAGE="adimpro")$eimg
  
  dim(yimg) <- dimg
  
  invisible(sqrt(ximg^2 + yimg^2))
}

shrink.image.old <- function(img, method = "gap", xt = img$dim[1], yt = img$dim[2], ratio = TRUE, compress=TRUE) {
  if(!check.adimpro(img)) {
    cat(" Consistency check for argument object failed (see warnings). object is returned.\n")
    return(invisible(img)) 
  }
  if(img$compressed) img <- decompress.image(img)
  if(img$type=="RAW") {
     warning("Can not shrink RAW images, return original")
     return(invisible(img))
  }
  type <- switch(img$type,rgb="rgb",yuv="csp",yiq="csp",hsi="csp",xyz="csp",greyscale="grey",RAW="grey")
  dimg <- img$dim

  #    check for discrepancies in xt, yt and ratio
  if(is.null(xt)||xt<10||xt>dimg[1]) xt <- dimg[1]
  if(is.null(yt)||yt<10||yt>dimg[2]) yt <- dimg[2]
  if(xt==dimg[1]&&yt==dimg[2]) {
     warning("specified target size does not define a shrinkage operation, return original")
     return(invisible(img))
  }
  if(ratio){
    xratio <- xt/dimg[1]
    yratio <- yt/dimg[2]
    if(xratio <= yratio) yt <- trunc(dimg[2]*xratio) else xt <- trunc(dimg[1]*yratio)
  }
  imethod <- switch(method,gap=1,mean=2,nearest=3,1)
  dv <- 3
  img$img <- switch(type,
                    rgb=array(.Fortran("shrnkrgb",
                                       as.integer(img$img),
                                       as.integer(dimg[1]),
                                       as.integer(dimg[2]),
                                       as.integer(dv),
                                       imgnew=integer(xt*yt*dv),
                                       as.integer(xt),
                                       as.integer(yt),
                                       integer(xt+1),
                                       integer(yt+1),
                                       as.integer(imethod),
                                       DUP=FALSE,
                                       PACKAGE="adimpro")$imgnew,c(xt,yt,dv)),
                    csp=array(.Fortran("shrnkcsp",
                                       as.double(img$img),
                                       as.integer(dimg[1]),
                                       as.integer(dimg[2]),
                                       as.integer(dv),
                                       imgnew=numeric(xt*yt*dv),
                                       as.integer(xt),
                                       as.integer(yt),
                                       integer(xt+1),
                                       integer(yt+1),
                                       as.integer(imethod),
                                       DUP=FALSE,
                                       PACKAGE="adimpro")$imgnew,c(xt,yt,dv)),
                    grey=matrix(.Fortran("shrnkgr",
                                         as.integer(img$img),
                                         as.integer(dimg[1]),
                                         as.integer(dimg[2]),
                                         imgnew=integer(xt*yt),
                                         as.integer(xt),
                                         as.integer(yt),
                                         integer(xt+1),
                                         integer(yt+1),
                                         as.integer(imethod),
                                         DUP=FALSE,
                                         PACKAGE="adimpro")$imgnew,xt,yt))
  img$dim <- c(xt,yt)
  if(!is.null(img$ni)) img$ni<-NULL
  if(!is.null(img$hmax)) img$hmax<-NULL
  if(!is.null(img$ni0)) img$ni0<-NULL

#  img$ni  has an incompatible dimension
  invisible(if(compress) compress.image(img) else img)
}

shrink.image <- function(img, method = "median", xt = img$dim[1], yt = img$dim[2], ratio = TRUE, compress=TRUE) {
  if(!check.adimpro(img)) {
    cat(" Consistency check for argument object failed (see warnings). object is returned.\n")
    return(invisible(img)) 
  }
  if(img$compressed) img <- decompress.image(img)
  if(img$type=="RAW") {
     warning("Can not shrink RAW images, return original")
     return(invisible(img))
  }
  type <- switch(img$type,rgb="color",yuv="color",yiq="color",hsi="color",xyz="color",greyscale="grey",RAW="grey")
  dimg <- img$dim

  #    check for discrepancies in xt, yt and ratio
  if(is.null(xt)||xt<10) xt <- dimg[1]
  if(is.null(yt)||yt<10) yt <- dimg[2]
 # if(xt==dimg[1]&&yt==dimg[2]) {
 #   warning("specified target size does not define a shrinkage operation, return original")
 #  return(invisible(img))
 #}
  if(ratio){
    xratio <- xt/dimg[1]
    yratio <- yt/dimg[2]
    if(xratio <= yratio) yt <- trunc(dimg[2]*xratio) else xt <- trunc(dimg[1]*yratio)
  }
  nz <- as.integer((dimg[1]/xt+2)*(dimg[2]/yt+2))
#  this is probably to much
  imethod <- switch(method,nearest=1,median=2,mean=3,1)
  dv <- 3
  img$img <- switch(type,
                    color=array(.Fortran("shrinkc",
                                       as.integer(img$img),
                                       as.integer(dimg[1]),
                                       as.integer(dimg[2]),
                                       imgnew=integer(xt*yt*dv),
                                       as.integer(xt),
                                       as.integer(yt),
                                       as.double(1.e-3),
                                       double(nz*dv),
                                       as.integer(nz),
                                       as.integer(imethod),
                                       DUP=FALSE,
                                       PACKAGE="adimpro")$imgnew,c(xt,yt,dv)),
                    grey=matrix(.Fortran("shrinkg",
                                       as.integer(img$img),
                                       as.integer(dimg[1]),
                                       as.integer(dimg[2]),
                                       imgnew=integer(xt*yt),
                                       as.integer(xt),
                                       as.integer(yt),
                                       as.double(1.e-3),
                                       double(nz),
                                       as.integer(nz),
                                       as.integer(imethod),
                                       DUP=FALSE,
                                       PACKAGE="adimpro")$imgnew,xt,yt))
  img$dim <- c(xt,yt)
  if(!is.null(img$ni)) img$ni<-NULL
  if(!is.null(img$hmax)) img$hmax<-NULL
  if(!is.null(img$ni0)) img$ni0<-NULL

#  img$ni  has an incompatible dimension
  invisible(if(compress) compress.image(img) else img)
}

rotate.image <- function(img,angle=90,compress=NULL) {
  # rotate image by 0, 90, 180 or 270 degrees
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  dimg <- img$dim
  if(is.null(compress)) compress <- img$compressed
  if(img$compressed) img <- decompress.image(img)
  if (img$type=="greyscale"||img$type=="RAW") {
    img$img <- switch(as.character(angle),
                      "0"=img$img,
                      "90"=t(img$img)[,dimg[1]:1],
                      "180"=img$img[dimg[1]:1,dimg[2]:1],
                      "270"=t(img$img)[dimg[2]:1,],
                      stop("Error: rotation by 0, 90, 180 or 270 degrees only"))
  } else {
    img$img <- switch(as.character(angle),
                      "0"=img$img,
                      "90"=aperm(img$img,c(2,1,3))[,dimg[1]:1,],
                      "180"=img$img[dimg[1]:1,dimg[2]:1,],
                      "270"=aperm(img$img,c(2,1,3))[dimg[2]:1,,],
                      stop("Error: rotation by 0, 90, 180 or 270 degrees only"))
  }
  if(!is.null(img$ni)) {
    img$ni <- switch(as.character(angle),
                      "0"=img$ni,
                      "90"=t(img$ni)[,dimg[1]:1],
                      "180"=img$ni[dimg[1]:1,dimg[2]:1],
                      "270"=t(img$ni)[dimg[2]:1,])
  }
  if(!is.null(img$xind)||!is.null(img$yind)) {
    xind <- img$xind
    yind <- img$yind
    img$xind <- NULL
    img$yind <- NULL
    if(!is.null(xind)) switch(as.character(angle),
                              "0"=img$xind<-xind,
                              "90"=img$yind<-xind[length(xind):1],
                              "180"=img$xind<-xind[length(xind):1],
                              "270"=img$yind<-xind) 
    if(!is.null(yind)) switch(as.character(angle),
                              "0"=img$yind<-yind,
                              "90"=img$xind<-yind,
                              "180"=img$yind<-yind[length(yind):1],
                              "270"=img$xind<-yind[length(yind):1]) 
    
  }
  img$dim <- dim(img$img)[1:2]
  img$rotate <- if(is.null(img$rotate)) angle/90 else (angle/90+img$rotate)%%4
  invisible(if(compress) compress.image(img) else img)
}

clip.image <- function(img,xind=NULL,yind=NULL,compress=NULL,...) {
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  #  if(is.null(compress))  compress <- img$compressed 
  #  if(img$compressed) img <- decompress.image(img)
  dimg <- img$dim
  ldimg <- switch(img$type,greyscale=2,RAW=2,3)
  if(!is.null(xind)||!is.null(yind)){
    if(is.numeric(xind)) xind <- as.integer(xind)
    if(is.null(xind)||!valid.index(xind,dimg[1])) xind <- 1:dimg[1] else xind <- (1:dimg[1])[xind]
    if(is.numeric(yind)) yind <- as.integer(yind)
    if(is.null(yind)||!valid.index(yind,dimg[2])) yind <- 1:dimg[2] else yind <- (1:dimg[2])[yind]
  } else {
    show.image(img,main="Identify the corners of the clipping region by left mouse clicks",...)
    cat("Identify the corners of the clipping region by left mouse clicks\n")
    coord <- locator(2,type="p")
    if(img$type=="RAW"){
#  require coordinates with same location of Bayer mask
    xa <- 2*(max(min(coord$x),1)%/%2)-1
    xe <- 2*(min(max(coord$x),dimg[1])%/%2)
    xind <- as.integer(xa:xe)
    ya <- 2*(max(min(coord$y),1)%/%2)-1
    ye <- 2*(min(max(coord$y),dimg[2])%/%2)
    yind <- as.integer(ya:ye)
    } else {
    xind <- as.integer(max(min(coord$x),1):min(max(coord$x),dimg[1]))
    yind <- as.integer(max(min(coord$y),1):min(max(coord$y),dimg[2]))
    }
    lines(c(min(xind),max(xind)),c(max(yind),max(yind)))
    lines(c(min(xind),max(xind)),c(min(yind),min(yind)))
    lines(c(min(xind),min(xind)),c(min(yind),max(yind)))
    lines(c(max(xind),max(xind)),c(min(yind),max(yind)))
  }
  if(img$compressed){
    dim(img$img) <- switch(img$type,greyscale=c(2,dimg),RAW=c(2,dimg),rgb=c(2,dimg,3),c(4,dimg,3))
    img$img <- as.vector(switch(img$type,greyscale=img$img[,xind,yind],RAW=img$img[,xind,yind],img$img[,xind,yind,]))
    if(!is.null(img$ni)) {
      dim(img$ni) <- c(4,dimg)
      img$ni <- as.vector(img$ni[,xind,yind])
    }
  } else {
    img$img <- switch(ldimg-1,img$img[xind,yind],img$img[xind,yind,])
    if(!is.null(img$ni)) img$ni <- img$ni[xind,yind]
  }
  img$dim <- c(length(xind),length(yind))
  if(is.null(img$xind)) img$xind <- xind else img$xind <- img$xind[xind]
  if(is.null(img$yind)) img$yind <- yind else img$yind <- img$yind[yind]
  show.image(img,main="Clipped image",...)
  invisible(img)
}


