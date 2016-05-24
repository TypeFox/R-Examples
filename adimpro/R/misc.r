.onLoad <- function(lib, pkg){
  # we need the path to Imagemagick
  mogrify <- Sys.which("mogrify")
  path2imagemagick <- strsplit(mogrify,"mogrify")[[1]][1]
  if(file.exists(mogrify) & (!is.na(path2imagemagick))) {
    Sys.setenv(ImageMagick=path2imagemagick)
  } else {
    warning("could not determine path to Imagemagick \n
    please set the correct path manually using \n
    'Sys.setenv(ImageMagick='path2imagemagick')'")
  }
  dcraw <- Sys.which("dcraw")
  if(!file.exists(dcraw)) packageStartupMessage("Reading RAW images requires to install dcraw, see \n
    http://cybercom.net/~dcoffin/dcraw/ for LINUX and http://www.insflug.org/raw/ 
    for MAC OS and Windows \n")
  adimpro.options()
}

adimpro.options <- function(xsize = NULL,
                            ysize = NULL){
  if (Sys.info()["sysname"] == "Linux") {
    ## we now have >which< on our system and test for >xdpyinfo< 
    test <- system("which xdpyinfo", ignore.stderr=TRUE, intern = TRUE)
    if (length(test) > 0) {
      resolution <- strsplit(system("xdpyinfo  | grep 'dimensions:'", intern=TRUE, ignore.stderr=TRUE), "x")
      if(length(resolution) > 0){
        resx <- strsplit(resolution[[1]][1]," ")[[1]]
        if (is.null(xsize)) xsize <- as.integer(resx[length(resx)])
        if (is.null(ysize)) ysize <- as.integer(strsplit(resolution[[1]][2]," ")[[1]][1])
      }
    }
  }
  if(is.null(xsize)||is.na(xsize)) xsize <- 1024
  if(is.null(ysize)||is.na(ysize)) ysize <- 768
  assign(".adimpro",list(xsize=xsize,ysize=ysize),pos=1)
  invisible(NULL)
}
#.adimpro <- adimpro.options()

check.adimpro <- function(object){
  # Returns true if object is of class adimpro, has
  # all necessary components and no contradicting information
  check <- NULL
  repeat{
    if(class(object)!="adimpro") {
      check <- 1
      break
    }
    if(is.null(object$type)) {
      check <- 4
      break
    } else type <- object$type
    if(!(type%in%c("rgb","hsi","yuv","yiq","xyz","greyscale","RAW"))) {
      check <- 4
      break
    }
    if(is.null(object$compressed)||!object$compressed){
      if(is.null(object$img)) {
        check <- 2
        break
      } else dimg <- dim(object$img)
      
      if(!(length(dimg)%in%c(2,3))) {
        check <- 2
        break
      }
      if(is.null(object$dim)) {
        check <- 3
        break
      } else dim <- dim(object$dim) 
      if(!all(dimg[1:2]==dim)) {
        check <- 3
        break
      }
      if(!(length(dimg)==switch(type,hsi=3,yuv=3,yiq=3,xyz=3,rgb=3,greyscale=2,RAW=2))) {
        check <- 5
        break
      }
      if(type %in% c("rgb","yuv","yiq","xyz","hsi") & (dimg[3]!=3)) {
        check <- 5
        break
      }
    }
    if(mode(object$gamma)!="logical") {
      check <- 6
      break
    } 
    if(is.null(object$depth)){
      check <- 7
      break
    } else depth <- object$depth
    if(!(depth %in% c("8bit","16bit"))){
      check <- 7
      break
    } 
    if(mode(object$wb)!="character") {
      check <- 8
      break
    } else wb <- object$wb
    if(!(wb %in% c("NONE","UNKNOWN","AUTO","CAMERA","USER","IMAGE","MAKE.IMAGE"))) {
      check <- 8
      break
    }
    if(wb == "USER" & is.null(object$whitep)) {
      check <- 9
      break
    }
    if(object$gammatype=="histogram"&&(is.null(object$hequal)||length(object$hequal)!=65536)) {
      check <- 10
      break
    }
    break
  }
  if(is.null(check)) adimpro <- TRUE else {
    adimpro <- FALSE
    warning(switch(check,
                   "Object is not of class adimpro",
                   "object$img is not a matrix or 3D-array",
                   "object$dim incompatible with dim(object$img)",
                   "object$type is not 'rgb', 'yuv', 'yiq', 'xyz', 'hsi', or 'greyvalue'",
                   "object$type does not correspond to dim(object$img)",
                   "object$gamma is not logical",
                   "object$depth is not '8bit' or '16bit'",
                   "object$wb is not one of 'NONE', 'Unknown', 'AUTO', 'CAMERA', 'Image', 'MAKE.IMAGE' or 'USER'",
                   "object$wb 'USER' but object$whitep not specified",
                   "object$gammatype is 'histogram' but object$hequal is invalid"))
  }
  adimpro
}

Spatialvar.gauss <- function(h,h0) {
  #   Calculates the factor of variance reduction obtained for Gaussian Kernel and bandwidth h in 
  #   case of colored noise that was produced by smoothing with Gaussian kernel and bandwidth h0
  #   Spatialvariance(lkern,h,h0,d)/Spatialvariance(lkern,h,1e-5,d) gives the 
  #   a factor for lambda to be used with bandwidth h 
  if(length(h)==1) h<-rep(h,2)
  ih<-trunc(4*h)
  ih<-pmax(1,ih)
  dx<-2*ih+1
  penl<-outer(dnorm(((-ih[1]):ih[1])/h[1]),dnorm(((-ih[2]):ih[2])/h[2]),"*")
  dim(penl)<-dx
  if(length(h0)==1) h0<-rep(h0,2)
  ih<-trunc(4*h0)
  ih<-pmax(1,ih)
  dx0<-2*ih+1
  x<- ((-ih[1]):ih[1])/h0[1]
  penl0<-outer(dnorm(((-ih[1]):ih[1])/h0[1]),dnorm(((-ih[2]):ih[2])/h0[2]),"*")
  dim(penl0)<-dx0
  penl0<-penl0/sum(penl0)
  dz<-dx+dx0-1
  z<-array(0,dz)
  for(i1 in 1:dx0[1]) for(i2 in 1:dx0[2]){
    ind1<-c(0:(i1-1),(dz[1]-dx0[1]+i1):dz[1]+1)
    ind1<-ind1[ind1<=dz[1]][-1]
    ind2<-c(0:(i2-1),(dz[2]-dx0[2]+i2):dz[2]+1)
    ind2<-ind2[ind2<=dz[2]][-1]
    z[-ind1,-ind2]<-z[-ind1,-ind2]+penl*penl0[i1,i2]
  }
  sum(z^2)/sum(z)^2
}

geth.gauss <- function(corr,step=1.01) {
  #   get the  bandwidth for lkern corresponding to a given correlation
  #  keep it simple result does not depend on d
  if (corr < 0.065) {
    h <- 1e-5
  } else { 
    h <- .3
    z <- 0
    while (z<corr) {
      h <- h*step
      z <- get.corr.gauss(h,interv=2)
    }
    h <- h/step
  }
  h
}

get.corr.gauss <- function(h,interv=1) {
  #   Calculates the correlation of 
  #   colored noise that was produced by smoothing with "gaussian" kernel and bandwidth h
  #   Result does not depend on d for "Gaussian" kernel !!
  h <- h*interv
  ih <- trunc(4*h+ 2*interv-1)
  dx <- 2*ih+1
  penl <- dnorm(((-ih):ih)/h)
  sum(penl[-(1:interv)]*penl[-((dx-interv+1):dx)])/sum(penl^2)
}

Varcor.gauss <- function(h) {
  #   Calculates a correction for the variance estimate obtained by (IQRdiff(y)/1.908)^2
  #   in case of colored noise that was produced by smoothing with lkern and bandwidth h
  h<-pmax(h,1e-5)
  if(h[1]<.25) vcg <- 1 else {
    if(length(h)<2) h<-rep(h,2)
    ih<-trunc(4*h)
    dx<-2*ih+1
    penl <- outer(dnorm(((-ih[1]):ih[1])/h[1]),dnorm(((-ih[2]):ih[2])/h[2]),"*")
    vcg <- 2*sum(penl)^2/sum(diff(penl)^2) 
  }
  vcg
}

valid.index <- function(ind,n) {
  #  Check if ind contains a valid subindex from (1:n)
  if(!is.numeric(ind)) return(FALSE)
  if(any(!is.integer(ind))) return(FALSE)
  if(prod(range(ind))<1) return(FALSE)
  if(max(abs(ind))>n) return(FALSE)
  TRUE
}

mask.create <- function(img,range1=c(0,1),range2=c(0,1),range3=c(0,1),locate=TRUE) {
  #  create a mask (to be used with awsimage)
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(img$compressed) img <- decompress.image(img)
  dimg <- img$dim
  if(locate){
    show.image(img)
    cat("select corners of rectangular region by mouse klick\n")
    coord <- locator(2)
    coord1 <- as.integer(pmin(pmax(coord$x,1),dimg[1]))
    coord2 <- as.integer(pmin(pmax(coord$y,1),dimg[2]))
  } else {
    coord1 <- c(1,1)
    coord2 <- dimg[1:2]
  }
  print(coord1)
  print(coord2)
  mask <- matrix(FALSE,dimg[1],dimg[2])
  xind <- min(coord1):max(coord1)
  yind <- min(coord2):max(coord2)
  mask[xind,yind] <- TRUE
  if(length(dimg)==2){
    mask[img$img<range1[1]|img$img>range1[2]] <- FALSE
  } else {
    mask[(img$img[,,1]<range1[1])|(img$img[,,1]>range1[2])] <- FALSE
    mask[(img$img[,,2]<range2[1])|(img$img[,,2]>range2[2])] <- FALSE
    mask[(img$img[,,3]<range3[1])|(img$img[,,3]>range3[2])] <- FALSE
  }
  mask
}

compress.image <- function(img){
  if(!check.adimpro(img)) {
    cat(" Consistency check for argument object failed (see warnings). object is returned.\"n")
    return(invisible(img)) 
  }
  if(is.null(img$compressed)||!img$compressed){
    type <- img$type
    dim(img$img) <- NULL
    size <- switch(type,rgb=2,greyscale=2,RAW=2,4)
    img$img <- writeBin(as.vector(img$img),raw(),size)
    if(!is.null(img$ni)) {
      dim(img$ni) <- NULL
      img$ni <- writeBin(as.vector(img$ni),raw(),4)
    }
    img$compressed <- TRUE
  }
  invisible(img)
}

decompress.image <- function(img){
  if(!check.adimpro(img)) {
    cat(" Consistency check for argument object failed (see warnings). object is returned.\"n")
    return(invisible(img)) 
  }
  if(!is.null(img$compressed)&&img$compressed){
    type <- img$type
    size <- switch(type,rgb=2,greyscale=2,RAW=2,4)
    what <- switch(type,rgb="integer",greyscale="integer",RAW="integer","numeric")
    nn <- prod(img$dim)*switch(type,greyscale=1,RAW=1,3)
    img$img <- readBin(img$img,what,nn,size,signed=FALSE)
    dim(img$img) <- switch(type,greyscale=img$dim,RAW=img$dim,c(img$dim,3))
    if(!is.null(img$ni)) {
      dim(img$ni) <- NULL
      img$ni <- readBin(img$ni,"numeric",prod(img$dim),4)
      dim(img$ni) <- img$dim
    }
    img$compressed <- FALSE
  }
  invisible(img)
}

getvofh2 <- function(bw,lkern){
.Fortran("getvofh2",
         as.double(bw),
         as.integer(lkern),
         vol=double(1),
         PACKAGE="adimpro")$vol
}
geth2 <- function(x,y,lkern,value,eps=1e-2){
.Fortran("geth2",
         as.double(x),
         as.double(y),
         as.integer(lkern),
         as.double(value),
         as.double(eps),
         bw=double(1),
         PACKAGE="adimpro")$bw
}

median1 <- function(x,tol=1e-8){
   if(!is.null(dim(x))&&dim(x)[2]==3) {
      z <- .Fortran("median3",
                    as.double(x),
                    as.integer(dim(x)[1]),
                    median=double(3),
                    as.double(tol),
                    PACKAGE="adimpro")$median
   } else {
      z <- .Fortran("median1",
                    as.double(x),
                    as.integer(length(x)),
                    median=double(1),
                    as.double(tol),
                    PACKAGE="adimpro")$median
   }
z
}

