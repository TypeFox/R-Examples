## TODO: the extraction of EXIF tags with identify does not work!

read.image <- function(filename, compress=TRUE) {
  convert.path <- paste(Sys.getenv("ImageMagick"),"convert",sep="")
  fileparts <- strsplit(filename,"\\.")[[1]]
  ext <- tolower(fileparts[length(fileparts)])
  name <- filename
  
  if (ext %in% c("tif","tiff","pgm","ppm","png","pnm","gif","jpg","jpeg","bmp")) {
    if (ext %in% c("tif","tiff","png","gif","jpg","jpeg","bmp")) {
      imgctype <- system(paste(Sys.getenv("ImageMagick"),"identify -format %r \"",filename,"\"",sep=""),intern=TRUE)
      if(regexpr("Gray",imgctype)!=-1){
      tmpfile <- paste(tempfile("pgm"),".pgm",sep="")
      tmpext <- "pgm"
      } else {
      tmpfile <- paste(tempfile("ppm"),".ppm",sep="")
      tmpext <- "ppm"
      }
      if (file.exists(filename)) {
            system(paste(convert.path," -compress None \"",filename,"\" \"",tmpfile,"\"",sep=""),wait=TRUE)
            filename <- tmpfile            
      } else {
        stop(paste("Error: file",filename,"does not exist!"))
      }
    }
    if(ext%in%c("pgm","ppm")) tmpext <- ext
    object <- list()
    
    if (tmpext == "pgm") {
      if (file.exists(filename)) {
        object$img <- read.pgm(filename)
      } else {
        stop("cannot find ",filename,"\n")
      }
      object$type <- "greyscale"
      object$cspace <- "greyscale"
    } else {
      if (file.exists(filename)) {
        object$img <- read.ppm(filename)
      } else {
        stop("cannot find ",filename,"\n")
      }
      object$type <- "rgb"
      object$cspace <- "sRGB"
    }
    object$depth <- attr(object$img, "depth")
    attr(object$img, "depth") <- NULL
    object$dim <- dim(object$img)[1:2]
    
    if (ext %in% c("tif","tiff","png","gif","jpg","jpeg")) file.remove(tmpfile)
    
    # in case of length(dim(img))==3  test if image contains same information in all channels
    if(length(dim(object$img)) == 3)
      if(sum(abs(range(object$img[,,1]-object$img[,,2]))+abs(range(object$img[,,1]-object$img[,,3])))==0) {
        object$img <- object$img[,,1]
        object$type <- "greyscale"
      }
    description <- paste(system(paste(Sys.getenv("ImageMagick"),"identify -format %c \"",name,"\" ",sep=""),intern=TRUE),collapse="\n",sep="")
    interpolation <- extract.info(description,"Interpolation")
    if(is.null(interpolation)) interpolation <- "unknown"
    object$interp <- interpolation
    gammatype <- extract.info(description,"Gammatype")
    if(is.null(gammatype)) gammatype <- "ITU"
    object$gammatype <- gammatype
    object$gamma <- gammatype %in% c("ITU","sRGB","CIE")
    whitep <- extract.info(description,"WhitePoint")
    if(is.null(whitep)) whitep <- "D65"
    object$whitep <- whitep
    wb <- extract.info(description,"WhiteBalance")
    if(is.null(wb)) wb <- "IMAGE"
    object$wb <- wb
    cspace <- extract.info(description,"Cspace")
    if(is.null(cspace)) cspace <- if(ext == "pgm") "greyscale" else "sRGB"
#  if cspace %in% c("xyz","hsi","yiq","yuv") we need to reconstruct the image
    if(cspace %in% c("xyz","hsi","yiq","yuv")){
       object$img <- object$img/65535
       object$type <- cspace
       compress <- FALSE
    }
    if(cspace == "yiq") {
       object$img[,,2] <- 1.2*object$img[,,2]-0.595716
       object$img[,,3] <- 1.06*object$img[,,3]-0.522591
    }
    if(cspace == "yuv") {
       object$img[,,2] <- object$img[,,2]-0.493
       object$img[,,3] <- 1.76*object$img[,,3]-0.877
    }
    object$cspace <- cspace 
    if(compress) {
      dim(object$img) <- NULL 
      object$img <- writeBin(as.integer(object$img),raw(),2)
    }
    file <- extract.info(description,"File")
    if(is.null(file)) file <- name
    object$file <- file
    object$xind <- extract.info(description,"xind")
    object$yind <- extract.info(description,"yind")
    ctype <- extract.info(description,"Type")
    if((!is.null(ctype))&&ctype=="RAW") {
       object$type <- "RAW"
       }
    object$description <- description
    object <- modify.info(object,add=FALSE)
# remove duplicated informmation from description
    object$compressed <- compress
    
    class(object) <- "adimpro"
    invisible(object)
    
  } else {
    warning(paste("Error: cannot handle file extension",ext))
    invisible(NULL)
  }
}

write.image <- function(img, file="tmp.ppm", max.x = NULL, max.y =NULL, depth=NULL,
       gammatype="ITU", whitep = NULL, temp = NULL, cspace = NULL, black=0, exposure=1) {
  
  convert.path <- paste(Sys.getenv("ImageMagick"),"convert",sep="")
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(is.null(cspace)) cspace <- img$cspace
  if(!(cspace %in% c("sRGB","Adobe","wGamut","kodak","xyz","yuv","yiq","hsi"))) {
       cspace <- "sRGB"
  }
  if(cspace %in% c("xyz","yuv","yiq","hsi")){
     gammatype <- "None"
  }
  if(img$compressed) img <- decompress.image(img)
  
  if(!is.null(max.x)||!is.null(max.y)) img <- shrink.image(img, xt=max.x, yt=max.y, method="nearest",compress=FALSE)
  
  fileparts <- strsplit(file,"\\.")[[1]]
  ext <- tolower(fileparts[length(fileparts)])
  
  # supported color depths
  if (is.null(depth)) {
    depth <- switch(img$depth,
                    "8bit" = 8,
                    "16bit" = 16,
                    8)
  }

  # image dimension
  dimg <- dim(img$img)
  
  # determine file name for PPM intermediate
  fileparts <- strsplit(file,"\\.")[[1]]
  ext <- tolower(fileparts[length(fileparts)])
  if (img$type == "greyscale") {
     tmpfile <- paste(tempfile("pgm"),".pgm",sep="")
  } else {
     tmpfile <- paste(tempfile("ppm"),".ppm",sep="")
  }
  
  # white balance
  dogamma <- !is.null(gammatype)&&gammatype!=img$gammatype
  docspace <- cspace!=img$cspace
  dowhite <- (!is.null(whitep)&&whitep!=img$whitep)||(!is.null(temp))
  doexposure <- !is.null(exposure)&&(exposure!=1||black!=0)
  if(dogamma||docspace||dowhite||doexposure){
  img <- adjust.image(img, gammatype=gammatype, cspace=cspace, 
                whitep=whitep, temp=temp, black=black, exposure=exposure, compress=FALSE) 
  }
  if(img$cspace %in% c("xyz","hsi")){
     img$img <- array(as.integer(65535*img$img),dim(img$img))
  }
  if(img$cspace=="yiq"){
     img$img[,,2] <- (img$img[,,2] + 0.595716)/1.2
     img$img[,,3] <- (img$img[,,3] + 0.522591)/1.06
     img$img <- array(as.integer(65535*img$img),dim(img$img))
  }
  if(img$cspace=="yuv"){
     img$img[,,2] <- (img$img[,,2] + 0.493)
     img$img[,,3] <- (img$img[,,3] + 0.877)/1.76
     img$img <- array(as.integer(65535*img$img),dim(img$img))
  }
  # rotate image appropriately
  pimg <- switch(img$type,
                 "greyscale" = img$img[,dimg[2]:1],
                 img$img[,dimg[2]:1,])
  
  # now write
  if(img$type != "unknown") {
    if (depth == 8) {
      pimg <- as.integer(pimg / 256)
      dim(pimg) <- dimg
      maximg <- 255
    } else {
      maximg <- 65535
    }
    ptype <- switch(img$type,
                    "greyscale" = "P5","P6")
    con <- file(tmpfile, "wb")
    writeChar(ptype,con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    writeChar(paste(dimg[1],dimg[2]),con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    writeChar(as.character(maximg),con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    close(con)
    
    if (img$type != "greyscale") pimg <- aperm(pimg,c(3,1,2))
    
    if (depth == 8) {
      con <- file(tmpfile, "ab")
      writeBin(as.vector(pimg), con, 1)
      close(con)      
    } else {
      con <- file(tmpfile, "ab")
      writeBin(as.vector(pimg), con, 2, endian="big")
      close(con)              
    }
    img <- modify.info(img,add=TRUE)
#  write additional information to image description
    if (tmpfile != file) {
          system(paste(convert.path,
           " -compress None -comment \"",img$description,"\" \"",
           tmpfile,"\" \"",file,"\"",sep=""),wait=TRUE)
          file.remove(tmpfile)
    } else {
          system(paste(Sys.getenv("ImageMagick"),
"mogrify -comment \"",img$description,"\" \"",tmpfile,"\"",sep=""),wait=TRUE)
    }
  } else {
    stop("Error: unknown colorspace type! exiting!")
  }
  invisible(NULL)
}
write.raw <- function(img, filename="tmp.png") {
  
  convert.path <- paste(Sys.getenv("ImageMagick"),"convert",sep="")
  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  fileparts <- strsplit(filename,"\\.")[[1]]
  if(!(tolower(fileparts[length(fileparts)]) %in% c("png"))){
     filename <- paste(filename,".png",sep="")
  }
  if(img$compressed) img <- decompress.image(img)
  depth <- 16
  dimg <- dim(img$img)
  
  # determine file name for PPM intermediate
    tmpfile <- paste(tempfile("pgm"),".pgm",sep="")
    pimg <- img$img[,dimg[2]:1]
  # now write
    maximg <- 65535
    ptype <- "P5"
    con <- file(tmpfile, "wb")
    writeChar(ptype,con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    writeChar(paste(dimg[1],dimg[2]),con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    writeChar(as.character(maximg),con,eos=NULL)
    writeBin(charToRaw("\n"),con)
    close(con)
    con <- file(tmpfile, "ab")
    writeBin(as.vector(pimg), con, 2, endian="big")
    close(con)
    img <- modify.info(img,add=TRUE)
          system(paste(convert.path,
       " -compress None -comment \"",img$description,"\" \"",
         tmpfile,"\" \"",filename,"\"",sep=""),wait=TRUE)
          file.remove(tmpfile)
  invisible(NULL)
}

show.image <- function (img, max.x = 1e+03, max.y =1e+03, gammatype = "ITU", whitep = NULL, temp = NULL,
                        cspace = "sRGB", black=0, exposure=1, channel=NULL, new = FALSE, ...) {

  if(!check.adimpro(img)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  
  if(cspace %in% c("grey","gray","grayscale")) cspace <- "greyscale" 
  dimg0 <- img$dim
  if(dimg0[1]>=max.x||dimg0[2]>=max.y){
  img <- shrink.image(img,xt=max.x,yt=max.y,method="nearest",compress=FALSE)
  } else if(img$compressed) img <- decompress.image(img)
  dimg <- img$dim
  
  
  # we need of course!
  di <- dim(img$img)
  
  # convert colorspace if neccessary
  if(!(cspace %in% c("sRGB","Adobe","wGamut","kodak"))) gammatype <- "None"
  if(toupper(img$type) == "RAW") {
#  avoid calling adjust.image 
     cspace <- img$cspace
     whitep <- temp <- NULL
     exposure <- 1
     black <- 0
     gammatype <- img$gammatype
  }
  if( cspace != img$cspace || !is.null(whitep) || !is.null(temp) ||exposure!=1||black!=0||gammatype!=img$gammatype) {
      img <- adjust.image(img,gammatype = gammatype, whitep = whitep, temp = temp,
                              cspace = cspace, black= black, exposure = exposure, compress=FALSE)
  }
#   if color.space %in% c("hsi","yuv","yiq","xyz") scale channels to provide maximum contrast
  if(img$type %in% c("yuv","yiq")){
     for(i in 1:3){
        img$img[,,i] <- as.integer(65535*(img$img[,,i]-min(img$img[,,i]))/(max(img$img[,,i])-min(img$img[,,i])))
     }
  }
  if(img$type %in% c("xyz","hsi")) img$img <- img$img*65535
   if(img$type=="hsi") img$img[,,1] <- (img$img[,,1]-min(img$img[,,1]))/diff(range(img$img[,,1]))*65535
   if (!is.null(channel)&&img$type!="greyscale") {
    if(channel %in% (1:3)){
      img$img <- img$img[,,channel]
      img$type <- "greyscale"
    }
  }
  if (new) X11()
    
  # now plot according to image type attribute
  switch(img$type,
         "greyscale" = show.greyscale(img,dimg0,...),
         "RAW" = show.greyscale(img,dimg0,...),
         "rgb" = show.rgb(img,dimg0,...),
         "hsi" = show.rgb(img,dimg0,...),
         "yuv" = show.rgb(img,dimg0,...),
         "yiq" = show.rgb(img,dimg0,...),
         "xyz" = show.rgb(img,dimg0,...),
         return())
  
}

show.rgb <- function(img,dimg0,...) {

  # check colorspace
  if(any(is.na(img$img))) cat("NA's in img$img\n")
  img$img[img$img < 0] <- as.integer(0)
  img$img[is.na(img$img)] <- as.integer(0)
  minimg <- min(img$img)
  maximg <- max(img$img)
  
  # should never be executed!!!
  if (maximg > 65535) {
   if(img$type=="rgb")    warning("Found values smaller than 0 or larger than 65535. Please check!")
#    img$img <- 65535*(img$img - minimg) / (maximg - minimg)
    img$img[img$img> 65535] <- 65535
  }
  
  # end check colorspace
  dimg <- img$dim
  
  # define an image z of same size as img and fill with increasing numbers
  z <- matrix(1:prod(dimg),nrow = dimg[1],ncol = dimg[2])
    
    
  # define the clor map such that for every pixel the color is set
  # according to the value of z
  #    color <-
  #      rgb(img$img[,,1]/65535,img$img[,,2]/65535,img$img[,,3]/65535)
  # we do _not_ use build-in R function rgb() due to memory usage!
  # define 0 to 255 in hexadecimal notation
  hex <- c(0:9, LETTERS[1:6])
  hex <- paste(hex[(0:255)%/%16+1],hex[(0:255)%%16+1],sep="")
  color <- paste("#",hex[img$img[,,1]%/%256+1],hex[img$img[,,2]%/%256+1],hex[img$img[,,3]%/%256+1],sep="")
    
  x <- seq(1,dimg0[1],length=dimg[1])
  y <- seq(1,dimg0[2],length=dimg[2])
  # display the image
  image(x, y, z, col = color, asp = 1, xlab="",ylab="", ...)
}

show.greyscale <- function(img,dimg0,...) {
  # if max.pixel is exceeded, shrink the image to a displayable size 
  dimg <- img$dim
  
  # check colorspace
  img$img[img$img < 0] <- 0
  img$img[is.na(img$img)] <- 0
  minimg <- min(img$img)
  maximg <- max(img$img)
  
  # should never be executed!!!
  if (maximg > 65535) {
     if(img$type=="rgb") warning("Found values smaller than 0 or larger than 65535. Please check!")
#    img$img <- (img$img - minimg) / (maximg - minimg)
     img$img[img$img> 65535] <- 65535
  }
  
  # end check colorspace
  # define an image z of same size as img and fill with increasing numbers
  z <- matrix(1:prod(dimg),nrow = dimg[1],ncol = dimg[2])


  color <- grey(img$img/65535)
    
  x <- seq(1,dimg0[1],length=dimg[1])
  y <- seq(1,dimg0[2],length=dimg[2])

  # display the image
  image(x,y,z,col=color, asp=1, xlab="",ylab="",...)
}

read.ppm <- function(filename) {
  # read header
  con <- file(filename,"r")
  type <- readLines(con,n=1)
  count <- 0
  while (TRUE) {
    nl <- readLines(con,n=1)
    if ((regexpr("^ *#", nl) != -1) || (regexpr("^ *$", nl) != -1)) {
      count <- count + 1
    } else {
      size <- nl
      break
    }
  }
  size <- as.integer(strsplit(size," ")[[1]])
  sizex <- as.integer(size[1])
  sizey <- as.integer(size[2])
  maxval <- as.integer(readLines(con,n=1))
  close(con)

  # read all
  if (type == "P3") {
    img <- as.integer(scan(filename,"int",sizex*sizey*3,skip=3+count,quiet=TRUE))
    if (maxval < 256) {
      img <- img * 256 # make "16bit"
      dd <- "8bit"
    } else {
      dd <- "16bit"      
    }
  } else if (type == "P6") {
    con <- file(filename,"rb")
    ttt <- readLines(con,n=3+count) # header again
    if (maxval > 255) {
      # endianess not clear. use quasi-standard
      img <- readBin(con,"int",n=sizex*sizey*3,2,signed=FALSE,endian="big") 
      dd <- "16bit"
    } else {
      img <- readBin(con,"int",n=sizex*sizey*3,1,signed=FALSE)
      img <- img * 256 # make "16bit"
      dd <- "8bit"
    }
    close(con)
  } else {
    stop("Cannot handle PPM type",type,"\n")
  }

  img <- as.integer(img)
  dim(img) <- c(3,sizex,sizey)
  img <- aperm(img,c(2,3,1))[,sizey:1,]

  attr(img,"depth") <- dd
  invisible(img)
}

read.pgm <- function(filename) {
  # read header
  con <- file(filename,"r")
  type <- readLines(con,n=1)
  count <- 0
  while (TRUE) {
    nl <- readLines(con,n=1)
    if ((regexpr("^ *#", nl) != -1) || (regexpr("^ *$", nl) != -1)) {
      count <- count + 1
    } else {
      size <- nl
      break
    }
  }
  size <- as.integer(strsplit(size," ")[[1]])
  sizex <- as.integer(size[1])
  sizey <- as.integer(size[2])
  maxval <- as.integer(readLines(con,n=1))
  close(con)

  # read all
  if (type == "P2") {
    img <- as.integer(scan(filename,"int",sizex*sizey,skip=3+count,quiet=TRUE))
    if (maxval < 256) {
      img <- img * 256 # make "16bit"
      dd <- "8bit"
    } else {
      dd <- "16bit"
    }
  } else if (type == "P5") {
    con <- file(filename,"rb")
    ttt <- readLines(con,n=3+count) # header again
    if (maxval > 255) {
      # endianess not clear. use quasi-standard
      img <- readBin(con,"int",n=sizex*sizey*3,2,signed=FALSE,endian="big") 
      dd <- "16bit"
    } else {
      img <- readBin(con,"int",n=sizex*sizey*3,1,signed=FALSE)
      img <- img * 256 # make "16bit"
      dd <- "8bit"
    }
    close(con)
  } else {
    stop("Cannot handle PGM type",type,"\n")
  }
  
  img <- as.integer(img)
  dim(img) <- c(sizex,sizey)
  img <- img[,sizey:1]
  
  attr(img,"depth") <- dd
  invisible(img)
}

read.raw <- function (filename,type="PPM",wb="CAMERA",cspace="Adobe",interp="Bilinear",maxrange=TRUE,rm.ppm=TRUE,compress=TRUE) {
# check the image to be a grey-value png that can be interpreted as
# containing RAW-data 
# otherwise just take it as a color image without gamma correction
    fileparts <- strsplit(filename,"\\.")[[1]]
    ext <- fileparts[length(fileparts)]
    if(tolower(ext) == "png"){
    img <- read.image(filename)
    if(img$type=="RAW") {
       img$depth <- "16bit"
       img$orientation <- "horizontal"
       img$interp <- "None"
       img$cspace <- cspace
       img$whitep <- "D65"
       img$wb <- "CAMERA"
       if(type!="RAW"){
          img <- develop.raw(img,"BILINEAR",maxrange=TRUE)
       }
    }
  return(invisible(img))
  }
  opt1 <- switch(toupper(type),PPM="-4",RAW="-4 -d",HALFSIZE="-h",INFO="-i -v","-4")
  opt2 <- if (opt1 == "-i -v") NULL else switch(toupper(wb),NONE=NULL,AUTO="-a",CAMERA="-w")
  opt3 <- switch(cspace,RAW="-o 0",sRGB=NULL,Adobe="-o 2",wGamut="-o 3",kodak="-o 4",XYZ="-o 5",NULL)
  opt4 <- switch(interp,Bilinear="-q 0",VNG="-q 2",AHD="-q 3",FourC="-f","-q 0")
  #  VNG seems to provide minimal spatial correlation
  tmpfile0 <- tempfile("raw")
  tmpfile <- paste(c(tmpfile0,ext),collapse=".")
  file.copy(filename,tmpfile)
  system(paste("dcraw", opt1, opt2, opt3, opt4, tmpfile))

  object <- list()
  
  if(opt1 != "-i -v") {
    if(opt1 == "-4 -d") {
      filenamep <- paste(tmpfile0,".pgm",sep="") 
      object$img <- read.pgm(filenamep)
      object$type <- "RAW"
      interp <- "None"
      cspace <- "CAMERA"
      wb <- "None"
    } else {
      filenamep <- paste(tmpfile0,".ppm",sep="") 
      if (file.exists(filenamep)) {
        object$img <- read.ppm(filenamep)
        } else {
        # ImageMagick under WINDOWS contains version of dcraw that
        # creates image.ppm instead of filename.ppm
        filename2 <- "image.ppm"
        if (file.exists(filename2)) {
          object$img <- read.ppm(filename2)
          filenamep <- filename2
        } else {
          stop("cannot find neither ",filename," nor ",filename2,"\n")
        }
      }
      object$type <- "rgb"
    }
    if(rm.ppm) file.remove(filenamep)
  }
  description <- system(paste("dcraw -i -v", filename),intern=TRUE)
  object$description  <- paste(description,collapse="\n",sep="")
  object$depth <- attr(object$img, "depth")
  attr(object$img, "depth") <- NULL
  object$dim <- dim(object$img)[1:2]
  if(object$type=="RAW") object$orientation <- if(object$dim[1]>=object$dim[2]) "horizontal" else "vertical"
  object$file <- filename
  object$interp <- interp
  object$cspace <- cspace
  object$whitep <- "D65"
  object$gamma <- FALSE
  object$gammatype <- "None"
  object$wb <- toupper(wb)
  if(compress) {
    dim(object$img) <- NULL 
    object$img <- writeBin(as.integer(object$img),raw(),2)
  }
  object$compressed <- compress
  class(object) <- "adimpro"
  invisible(object)
}

  
plot.adimpro <- function(x, new = FALSE, gammatype=NULL, cspace=NULL, whitep=NULL, temp=NULL, black=0, exposure=1,...) {
  if(!check.adimpro(x)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if(x$compressed) x <- decompress.image(x)
  if(!(is.null(gammatype)&&is.null(cspace)&&is.null(whitep)&&is.null(temp)&&exposure==1&&black==0)) {
     x <- adjust.image(x, gammatype=gammatype, cspace=cspace, whitep=whitep, temp=temp, black=black,exposure=exposure,compress=FALSE)
  }
  if ("ni" %in% names(x)) {
    aws <- TRUE
    ni <- make.image(65535* x$ni / x$ni0, compress=FALSE, gammatype="None")
  } else {
    aws <- FALSE
  }
  if (x$type == "rgb") {
    col1 <- "red"
    col2 <- "green"
    col3 <- "blue"
    tt <- c("Red channel","Green channel","Blue channel")
    xlim1 <- c(0,65535)
    xlim2 <- c(0,65535)
    xlim3 <- c(0,65535)
  } else if (x$type == "hsi") {
    col1 <- hsv(0:99/99,rep(1,100),rep(1,100))
    col2 <- hsv(rep(1,100),0:99/99,rep(1,100))
    col3 <- grey(0:99/99)
    tt <- c("Hue","Saturation","Intensity")
    xlim1 <- c(0,1)
    xlim2 <- c(0,1)
    xlim3 <- c(0,1)
  } else if (x$type == "yuv") {
    col1 <- grey(0:99/99)
    col2 <- "green1"
    col3 <- "green4"
    tt <- c("Intensity","U channel","V channel")
    xlim1 <- c(0,1)
    xlim2 <- c(-.436,.436)
    xlim3 <- c(-0.615,0.615)
  } else if (x$type == "yiq") {
    col1 <- grey(0:99/99)
    col2 <- "blue1"
    col3 <- "blue4"
    tt <- c("Intensity","I channel","Q channel")
    xlim1 <- c(0,1)
    xlim2 <- c(-.596,.596)
    xlim3 <- c(-.523,.523)
  } else if (x$type == "xyz") {
    col1 <- grey(0:99/99)
    col2 <- "blue1"
    col3 <- "blue4"
    tt <- c("X channel","Y channel","Z channel")
    xlim1 <- range(x$img[,,1])
    xlim2 <- range(x$img[,,2])
    xlim3 <- range(x$img[,,3])
  } else if (x$type == "greyscale") {
    col <- grey(0:99/99)
    tt <- "Grey level"
    xlim <- c(0,65535)
  } else {
    stop("Really dont know what to do with this image type! Sorry!")
  }
  
  if (x$type == "greyscale") {
    if (aws) { 
      if (new) X11(width=5,height=5)
      oldpar <- par(mfrow = c(2,2),mar=c(3,3,3,1))
      on.exit(par(oldpar))
    } else {
      if (new) X11(width=7,height=3)
      oldpar <- par(mfrow = c(1,3),mar=c(3,3,3,1))
      on.exit(par(oldpar))
    }
    hist(x$img,100,col=col,border=col,xlim=xlim,main=tt)
    show.image(x,max.x=400,max.y=400,gammatype=if(is.null(gammatype)) "ITU" else gammatype)
    
    plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",frame.plot=FALSE)
    text(0,1,paste("Image file",x$file),pos=4)
    text(0,0.9,paste("Min",min(x$img)),pos=4)
    text(0,0.8,paste("Max",max(x$img)),pos=4)
    
    if (aws) show.image(ni,max.x=400,max.y=400)
    
  } else {
    if (new) X11(width=7,height=5)
    oldpar <- par(mfrow = c(2,3),mar=c(3,3,3,1))
    on.exit(par(oldpar))
    
    hist(x$img[,,1],100,col=col1,border=col1,xlim=xlim1,main=tt[1])
    hist(x$img[,,2],100,col=col2,border=col2,xlim=xlim2,main=tt[2])
    hist(x$img[,,3],100,col=col3,border=col3,xlim=xlim3,main=tt[3])
    show.image(x,max.x=400,max.y=400, gammatype=if(is.null(gammatype)) "ITU" else gammatype)
    
    plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",frame.plot=FALSE)
    text(0,1,paste("Image file",x$file),pos=4)
    text(0,0.9,paste("Min",min(x$img)),pos=4)
    text(0,0.8,paste("Max",max(x$img)),pos=4)
    text(0,0.7,paste("Colorspace",x$cspace),pos=4)
    if(is.character(x$whitep)) wp <- x$whitep else wp <- paste("Temp=",temp)
    text(0,0.6,paste("White point",wp, paste("(",paste(signif(whitepoint(x$whitep),4),collapse=","),")")),pos=4)
    if(black!=0) text(0,0.5,paste("Black",black),pos=4)
    if(exposure!=1) text(0,0.4,paste("Exposure",exposure),pos=4)
    text(0,0.3,paste("Gamma",x$gammatype),pos=4)
    
    if (aws) show.image(ni,max.x=400,max.y=400)
  }

}

summary.adimpro <- function(object, ...) {
  if(!check.adimpro(object)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  #  if(object$compressed) object <- decompress.image(object)
  cat("         Image file:", object$file,"\n") 
  if(!is.null(object$xind))
    cat("horizontal clipping:", min(object$xind),":",max(object$xind),"\n")
  if(!is.null(object$yind))
    cat("  vertical clipping:", min(object$yind),":",max(object$yind),"\n")
  if(!is.null(object$rotate))
    cat(" Image rotated by:", object$rotate*90,"degrees\n")
  cat("    Image dimension:", object$dim,"\n")
  cat("        Color space:", object$type)
  if(object$type=="rgb") cat("(", object$cspace, ")")
  cat("\n")
  cat("        Color depth:", object$depth,"\n")
  cat("   Gamma correction:", object$gamma," Type:", object$gammatype,"\n")
  cat("        White point:", object$whitep,"\n")
  if(!object$compressed) cat("              Range:", as.integer(range(object$img)),"\n")
  if(object$compressed)  cat("   Compressed image\n")
  if (!is.null(object$hmax))
    cat("     max. bandwidth:", signif(object$hmax,2),"\n")
  if (!object$compressed && !is.null(object$ni))
    cat("mean rel adaptation:", signif(mean(object$ni)/object$ni0,2),"\n")
  if (!is.null(object$scorr))
    cat("spatial correlation:", object$scorr,"\n")
  if (!is.null(object$chcorr))
    cat("channel correlation:", object$chcorr,"\n")
  if (!is.null(object$varcoef))
    cat("est. variance param:", object$varcoef,"\n")
  if(!is.null(object$description)) {
    cat("\nEXIF-Information:\n")
    cat(object$description,"\n")
  }
}

make.image <- function(x, compress=TRUE, gammatype="None", whitep="D65", cspace="Adobe", scale="Original", xmode="RGB"){
  if(gammatype %in% c("ITU","sRGB","CIE")) {
      gamma <- TRUE
  } else {
      gamma <- FALSE
      gammatype <- "None"
  }
  dimg <- dim(x)
  if(is.null(dimg) || !(length(dimg) %in% 2:3)) return(warning("x is not an array of appropriate dimensions."))
  if(diff(range(x))==0) {
     xmode=="RGB"
  }  
  if(xmode=="RGB"){
  if(diff(range(x))==0) {
      x <- 0*x 
  } else {
  if(scale=="Maxcontrast"){
     if(length(dimg)==2){
        x <- 65535*(x-min(x))/(max(x)-min(x))
     } else {
        for(i in 1:dimg[3]) x[,,i] <- 65535*(x[,,i]-min(x[,,i]))/(max(x[,,i])-min(x[,,i]))
     }
  } else {
     if(min(x) < 0) x <- (x-min(x))/(max(x)-min(x))*max(x)
     if( max(x) <= 1) x <- 65535 * x
     if( max(x) > 65535 ) x <- x/max(x)*65535
  } 
  }
  dim(x) <- NULL
  x <- if(compress) writeBin(as.integer(x),raw(),2) else array(as.integer(x),dimg)
  img <- list(img=x, type=switch(length(dimg)-1,"greyscale","rgb"),depth="16bit",
              dim=dimg[1:2], gamma=gamma, gammatype=gammatype,whitep=whitep,description="",
              cspace=switch(length(dimg)-1,"greyscale",cspace), file="artificial",wb="MAKE.IMAGE",compressed=compress)
  class(img) <- "adimpro"
  } else {
# interprete components in x as HSI
     if(min(x[,,1])<0||max(x[,,1])>1) {
        x[,,1] <- (x[,,1]-min(x[,,1]))/(max(x[,,1])-min(x[,,1]))
     }
     if(min(x[,,2])<0||max(x[,,2])>1) {
        x[,,2] <- (x[,,2]-min(x[,,2]))/(max(x[,,2])-min(x[,,2]))
     }
     if(min(x[,,3])<0||max(x[,,3])>1) {
        x[,,3] <- (x[,,3]-min(x[,,3]))/(max(x[,,3])-min(x[,,3]))
     }
  img <- list(img=x, type="hsi",dim=dimg[1:2],depth="16bit", gamma=gamma, wb="UNKNOWN", file="artificial",compressed=FALSE, gammatype=gammatype, whitep=whitep,description="")
  class(img) <- "adimpro"
  img <- hsi2rgb(img)
  }
  invisible(img) 
}

extract.ni <- function (object, gammatype = "ITU", compress=TRUE) {
  if (!check.adimpro(object)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (object$compressed)
    object <- decompress.image(object)
  if ("ni" %in% names(object)) {
    aws <- TRUE
    ni <- make.image(65535 * object$ni/object$ni0, compress = compress, gammatype = gammatype)
  }
  else {
    stop("image was not processed by awsimage or awspimage")
  }
  invisible(ni)
}

extract.image <- function (object) {
  if (!check.adimpro(object)) {
    stop(" Consistency check for argument object failed (see warnings).\n")
  }
  if (object$compressed)
    object <- decompress.image(object)
  invisible(object$img)
}

## TODO: This function should be generalized to arbitrary Tags and should still work for RAW and JPG

extract.info <- function(object, what="Bayer"){
  if (class(object)=="adimpro") {
      description <- object$description
  } else {
description <- as.character(object)
}
if(!is.null(description)){
erg<-switch(what,
            Bayer=substr(strsplit(description,"Filter pattern: ")[[1]][2],1,4),
            Daymulti=strsplit(description,"Daylight multipliers: ")[[1]][2],
            Cammulti=strsplit(description,"Camera multipliers: ")[[1]][2],
            Camera=strsplit(description,"Camera: ")[[1]][2],
            Isize=strsplit(description,"Image size:  ")[[1]][2],
            Osize=strsplit(description,"Output size: ")[[1]][2],
            File=strsplit(description,"File: ")[[1]][2],
            Interpolation=strsplit(description,"Interpolation: ")[[1]][2],
            Gammatype=strsplit(description,"Gammatype: ")[[1]][2],
            WhiteBalance=strsplit(description,"White Balance: ")[[1]][2],
            WhitePoint=strsplit(description,"White Point: ")[[1]][2],
            Cspace=strsplit(description,"Color Space: ")[[1]][2],
            Type=strsplit(description,"Type: ")[[1]][2],
            xind=strsplit(description,"Xind: ")[[1]][2],
            yind=strsplit(description,"Yind: ")[[1]][2],
            NULL
           )
if(!is.null(erg)&&is.na(erg)) erg <- NULL
# case of no information
if(!is.null(erg)) {
   erg <- strsplit(erg,"\n")[[1]][1]
   if(what%in%c("xind","yind")) {
      erg <- strsplit(erg,":")[[1]]
      erg <- as.integer(erg[1]):as.integer(erg[2])
   }
   if(what %in% c("Daymulti","Cammulti")){
      zz <- textConnection(erg)
      erg <- scan(zz, numeric(1))
      close(zz)
   }
}
} else {
erg <- NULL
}
erg
}
modify.info <- function(img,add=TRUE){
if(!is.null(img$description)){
z <- textConnection(img$description,"r")
description <- readLines(z)
close(z)
if(add){
description <- c(description,paste("File:",img$file),
                             paste("Interpolation:",img$interp),
                             paste("Gammatype:",img$gammatype),
                             paste("White Balance:",img$wb),
                             paste("White Point:",img$whitep),
                             paste("Type:",img$type),
                             paste("Color Space:",img$cspace),
                             if(!is.null(img$xind)) paste("Xind:",paste(min(img$xind),max(img$xind),sep=":")),
                             if(!is.null(img$yind)) paste("Yind:",paste(min(img$yind),max(img$yind),sep=":")))
} else {
#  remove components that are stored elsewhere
description <- strsplit(img$description,"File: ")[[1]][1]
}
}
img$description <- paste(description,collapse="\n",sep="")
img
}



develop.raw <- function(object,method="BILINEAR",wb=c(1,1,1),maxrange=TRUE,compress=TRUE){
#
#   converts Sensor data into RGB-images 
#   white balance is applied as a correction on the sensor data
#   in contrast to other functions where white balance is done in XYZ
#
  method <- toupper(method)
  if(!(method %in% c("FULL","HALF","BILINEAR","MEDIAN16","MEDIAN4"))) stop("Method not implemented")
  if(object$type!="RAW") stop("object does not contain RAW sensor data, 
                    please read the image by read.raw(filename,type=''RAW'')")
  if(object$compressed) object <- decompress.image(object)
  dimg <- dim(object$img)
  n1 <- dimg[1]
  n2 <- dimg[2]
  bayer <- switch(extract.info(object),RGGB=1,GRBG=2,BGGR=3,GBRG=4)
  if(extract.info(object,"Isize")!=extract.info(object,"Osize")) bayer <- bayer+1
  bayer <- (bayer-1)%%4+1
  if(!is.null(object$rotate)) {
      bayer <- object$rotate+bayer 
      bayer <- (bayer-1)%%4+1
  }
  if(any(wb!=1)){
# White balance if specified
     object$img <- matrix(.Fortran("wbalance",
                              sensor=as.integer(object$img),
                              as.integer(n1),
                              as.integer(n2),
                              as.double(wb),
                              as.integer(bayer),
                              DUP=FALSE,
                              PACKAGE="adimpro")$sensor,n1,n2)
  }
  if(maxrange){
     minimg <- min(object$img)
     rangeimg <- max(object$img)-minimg
     object$img <- matrix(as.integer((object$img-minimg)/rangeimg*65535),n1,n2)
  }
  h1 <- switch(method,FULL=n1-4,HALF=n1%/%2-1,MEDIAN16=n1-6,MEDIAN16=n1-2,n1)
  h2 <- switch(method,FULL=n2-4,HALF=n2%/%2-1,MEDIAN16=n2-6,MEDIAN16=n1-2,n2)
  theta <- array(switch(method,
                   FULL=.Fortran("fullsize",
                   as.integer(object$img),
                   theta=integer(h1*h2*3),
                   as.integer(n1),
                   as.integer(n2),
                   as.integer(h1),
                   as.integer(h2),
                   as.integer(bayer),
                   DUP=FALSE,
                   PACKAGE="adimpro")$theta,
                   HALF=.Fortran("halfsize",
                   as.integer(object$img),
                   theta=integer(h1*h2*3),
                   as.integer(n1),
                   as.integer(n2),
                   as.integer(h1),
                   as.integer(h2),
                   as.integer(bayer),
                   DUP=FALSE,
                   PACKAGE="adimpro")$theta,
                   BILINEAR=.Fortran("indemos4",
                   as.integer(object$img),
                   theta=integer(n1*n2*3),
                   as.integer(n1),
                   as.integer(n2),
                   as.integer(bayer),
                   as.integer(rep(1,3*n1*n2)),
                   integer(3*n1*n2),
                   DUP=FALSE,
                   PACKAGE="adimpro")$theta,
                   MEDIAN16=.Fortran("demmed16",
                   as.integer(object$img),
                   theta=integer(h1*h2*3),
                   as.integer(n1),
                   as.integer(n2),
                   as.integer(h1),
                   as.integer(h2),
                   as.integer(bayer),
                   DUP=FALSE,
                   PACKAGE="adimpro")$theta,
                   MEDIAN4=.Fortran("demmed4",
                   as.integer(object$img),
                   theta=integer(h1*h2*3),
                   as.integer(n1),
                   as.integer(n2),
                   as.integer(h1),
                   as.integer(h2),
                   as.integer(bayer),
                   DUP=FALSE,
                   PACKAGE="adimpro")$theta),c(h1,h2,3))
  n <- h1*h2
  out.cam <- cam2rgbmat(object)
  object$img <- array(.Fortran("cam2rgb",
                                as.integer(theta),
                                as.integer(n),
                                as.double(out.cam),
                                theta=integer(n*3),
                                DUP=FALSE,
                                PACKAGE="adimpro")$theta,c(h1,h2,3))
  object$type <- "rgb"
  object$dim <- c(h1,h2)
  object$interp <- method
  if(method=="HALF"){
     if(!is.null(object$xind)) object$xind<-unique(object$xind[-1]%/%2)
     if(!is.null(object$yind)) object$yind<-unique(object$yind[-1]%/%2)
  }
  invisible(if(compress) compress.image(object) else object)
}

adimpro2biOps <- function(img) {
  if (!("adimpro" %in% class(img))) stop("No adimpro object!")
  
  if (img$compressed) {
    nimg <- extract.image(img)
  } else {
    nimg <- img$img
  }
  if (length(dim(nimg)) == 3) { # RGB image
    nimg <- aperm(nimg, c(2,1,3))
    type <- "rgb"
   nimg <- nimg[dim(nimg)[1]:1,,]
  } else if (length(dim(nimg)) == 2) { # greyscale
    nimg <- t(nimg)
    type <- "grey"
    nimg <- nimg[dim(nimg)[1]:1,]
  } else {
    stop("ERROR: Image data has dimension",dim(nimg),"do not know what to do! Expecting (x,y,3) or (x,y)!")
  }
  nimg <- nimg/256
  
  attr(nimg, "type") <- type
  class(nimg) <- c("imagedata", class(nimg))
  invisible(nimg)
}

biOps2adimpro <- function(img,
                          gammatype = "ITU",
                          whitep    = "D65",
                          cspace    = "sRGB") {
  if (!("imagedata" %in% class(img))) stop("No imagedata (biOps) object!")

  if (length(dim(img)) == 3) { # RGB image
    img <- aperm(img, c(2,1,3))
    type <- "rgb"
    img <- img[,dim(img)[2]:1,]
  } else if (length(dim(img)) == 2) { # greyscale
    img <- t(img)
    type <- "greyscale"
    img <- img[,dim(img)[2]:1]
  } else {
    stop("ERROR: Image data has dimension",dim(img),"do not know what to do! Expecting (x,y,3) or (x,y)!")
  }
  img <- img/255 #  biOps seems to expect range [0,255]
  
  # don't know which default values!!
  invisible(make.image(img,
                       compress  = TRUE,
                       gammatype = gammatype,
                       whitep    = whitep,
                       cspace    = cspace,
                       scale     = "Original",
                       xmode     = "RGB"))
}
