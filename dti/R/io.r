################################################################
#                                                              #
# Section for data reading functions                           #
#                                                              #
################################################################

dtiData <- function(gradient,imagefile,ddim,bvalue=NULL,xind=NULL,yind=NULL,zind=NULL,level=0,mins0value=1,maxvalue=32000,voxelext=c(1,1,1),orientation=c(0L,2L,5L),rotation=diag(3)) {
  args <- list(sys.call())
  if (any(sort((orientation)%/%2) != 0:2)) stop("invalid orientation \n")
  if (dim(gradient)[2]==3) gradient <- t(gradient)
  if (dim(gradient)[1]!=3) stop("Not a valid gradient matrix")
  ngrad <- dim(gradient)[2]
  s0ind <- (1:ngrad)[apply(abs(gradient),2,max)==0] 
  if(length(s0ind)>0){
     d <- diag(t(gradient[,-s0ind])%*%gradient[,-s0ind])
     gradient[,-s0ind] <- t(t(gradient[,-s0ind])/sqrt(d))
     if (is.null(bvalue)){
        bvalue <- rep(1,ngrad)
        bvalue[s0ind] <- 0
     } else {
       if(length(bvalue)!=ngrad || max(bvalue[s0ind]) > 10*min(bvalue[-s0ind])) 
         stop("invalid b-values")
     }   
  } else {
## some datasets have small b-values assigned to unweighted data, e.g. HCP
## define s0ind as set of minimal b-values
     if(is.null(bvalue)) stop("need bvalues")
     if(length(bvalue)!=ngrad) stop("invalid b-values")
     s0ind <- (1:ngrad)[bvalue==min(bvalue)]
     cat("identify S0 images with b-value",min(bvalue),"\n")
     if(length(s0ind)==ngrad) stop("nonidentifability, all b-values are equal")
     d <- diag(t(gradient)%*%gradient)
     gradient <- t(t(gradient)/sqrt(d))
  }
  if (!(file.exists(imagefile))) stop("Image file does not exist")
  cat("Start Data reading",format(Sys.time()), "\n")
  zz <- file(imagefile,"rb")
  #  si now contains all images (S_0 and S_I), ngrad includes 
  #  number of zero gradients
  
  if (is.null(xind)) xind <- 1:ddim[1]
  if (is.null(yind)) yind <- 1:ddim[2]
  if (is.null(zind)) zind <- 1:ddim[3]
  si <- numeric()
  for (grad in 1:ngrad) {
    sitemp <- readBin(zz,"integer",prod(ddim),2,FALSE)
    dim(sitemp) <- ddim
    si <- c(si,sitemp[xind,yind,zind])
    cat(".")
  }
  close(zz)
  dim(si) <- c(length(xind),length(yind),length(zind),ngrad)
  dimsi <- dim(si)
  
  cat("Data successfully read",format(Sys.time()), "\n")
  
  #
  #   set correct orientation
  #
  xyz <- (orientation)%/%2+1
  swap <- orientation%%2
  if(any(xyz!=1:3)) {
    abc <- 1:3
    abc[xyz] <- abc
    si <- aperm(si,c(abc,4))
    swap[xyz] <- swap
    voxelext[xyz] <- voxelext
    dimsi[xyz] <- dimsi[1:3]
    ddim[xyz] <- ddim[1:3]
    gradient[xyz,] <- gradient
  }
  if(swap[1]==1) {
    si <- si[dimsi[1]:1,,,] 
    gradient[1,] <- -gradient[1,]
  }
  if(swap[2]==1) {
    si <- si[,dimsi[2]:1,,]  
    gradient[2,] <- -gradient[2,]
  }
  if(swap[3]==0) {
    si <- si[,,dimsi[3]:1,]    
    gradient[3,] <- -gradient[3,]
  }
  #
  #   orientation set to radiological convention
  #
  si <- .Fortran("initdata",
                 si=as.double(si),
                 as.integer(dimsi[1]),
                 as.integer(dimsi[2]),
                 as.integer(dimsi[3]),
                 as.integer(dimsi[4]),
                 as.double(maxvalue),
                 PACKAGE="dti")$si
  #  this replaces the content off all voxel with elements <=0 or >maxvalue by 0
  if(all(si==as.integer(si))) si <- as.integer(si)
  ##  reduce memory requirements if possible
  dim(si) <- dimsi
  level <- max(mins0value,level*mean(si[,,,s0ind][si[,,,s0ind]>0])) # set level to level*mean  of positive s_0 values
  ddim0 <- as.integer(ddim)
  ddim <- as.integer(dim(si)[1:3])
  
  cat("Create auxiliary statistics",format(Sys.time()), " \n")
  rind <- replind(gradient)
  design <- create.designmatrix.dti(gradient)
  invisible(new("dtiData",
                call = args,
                si     = si,
                gradient = gradient,
                bvalue = bvalue,
                btb    = sweep( design, 2, bvalue, "*"),
                ngrad  = ngrad, # = dim(btb)[2]
                s0ind  = s0ind, # indices of S_0 images
                replind = rind,
                ddim   = ddim,
                ddim0  = ddim0,
                xind   = xind,
                yind   = yind,
                zind   = zind,
                level  = level,
                sdcoef = c(1,0,level,quantile(si[,,,s0ind],.98)),
                voxelext = voxelext,
                orientation = as.integer(c(0,2,5)), #   orientation set to radiological convention
                rotation = rotation,
                source = imagefile)
  )
}

############


## TODO: what about siemens mosaic?
##       the loop over files contains many not neccessary operations (voxelext etc)
readDWIdata <- function(gradient, dirlist, 
                        format = c("DICOM", "NIFTI", "ANALYZE", "AFNI"), 
                        nslice = NULL, order = NULL, bvalue = NULL,
                        xind = NULL, yind = NULL, zind = NULL,
                        level = 0, mins0value = 1, maxvalue = 32000,
                        voxelext = NULL, orientation = c(0L, 2L, 5L), rotation = NULL,
                        pattern = NULL,
                        SPM2 = TRUE,
                        verbose = FALSE) {
  
  args <- list(sys.call())
  
  ## basic consistency checks
  format <- match.arg(format)
  if ((format == "DICOM") & is.null(nslice))
    stop("readDWIdata: Cannot handle DICOM folders without specifying number of slices nslice!")
  if (length(dim(gradient)) != 2) stop("readDWIdata: Not a valid gradient matrix.")
  if (dim(gradient)[2] == 3) gradient <- t(gradient)
  if (dim(gradient)[1] != 3) stop("readDWIdata: Not a valid gradient matrix.")
  ngrad <- dim(gradient)[2]
  s0ind <- (1:ngrad)[apply(abs(gradient), 2, max) == 0]
  if (length(s0ind) > 0) {
     d <- diag(t(gradient[,-s0ind])%*%gradient[,-s0ind])
     gradient[,-s0ind] <- t(t(gradient[,-s0ind])/sqrt(d))
     if (is.null(bvalue)){
        bvalue <- rep(1,ngrad)
        bvalue[s0ind] <- 0
     } else {
        if(length(bvalue)!=ngrad || max(bvalue[s0ind]) > 10*min(bvalue[-s0ind])) 
           stop("invalid b-values")
     } 
  } else {
## some datasets have small b-values assigned to unweighted data, e.g. HCP
## define s0ind as set of minimal b-values
     if(is.null(bvalue)) stop("need bvalues")
     if(length(bvalue)!=ngrad) stop("invalid b-values")
     s0ind <- (1:ngrad)[bvalue==min(bvalue)]
     cat("identify S0 images with b-value",min(bvalue),"\n")
     if(length(s0ind)==ngrad) stop("nonidentifability, all b-values are equal")
     d <- diag(t(gradient)%*%gradient)
     gradient <- t(t(gradient)/sqrt(d))
  }
  if(length(dirlist)==1&!is.na(isdir <- file.info(dirlist)$isdir)&!isdir){
    ## dirlist contains a filename rather than a list of directories
    filelist <- dirlist
    dirlist <- NULL
  } else {
    ## generate file list in specified order
    filelist <- NULL
    for (dd in dirlist) filelist <- c(filelist, list.files(dd, full.names = TRUE, pattern = pattern))
  }
  if ( length( filelist) == 0) stop( "readDWIdata: empty directories or directories do not exist!")
  if (format == "DICOM") {
    if (is.null(zind)) zind <- 1:nslice
    if (length(filelist) != ngrad * nslice)
      stop("readDWIdata: Number of found files (", length(filelist), ") does not match ngrad*nslice")
    if (is.null(order)) {
      order <- 1:(ngrad * nslice)
    } else {
      if (length(order) != ngrad * nslice)
        stop("readDWIdata: Length of order vector does not match ngrad*nslice")
    }
    dim(order) <- c(nslice, ngrad)
    order <- order[zind, ]
    dim(order) <- NULL
    filelist <- filelist[order]
  } else {
    if (format == "ANALYZE") filelist <- unlist(strsplit(filelist[regexpr("\\.hdr$", filelist) != -1], "\\.hdr"))
    if (format == "AFNI") filelist <- filelist[regexpr("\\.HEAD$", filelist) != -1]
    if (format == "NIFTI") {
      if ((length(filelist) == 2 * ngrad)) filelist <- unlist(strsplit(filelist[regexpr("\\.hdr$", filelist) != -1], "\\.hdr"))
      if ((length(filelist) != ngrad) & (length(filelist) != 1))
        stop("readDWIdata: Number of files (", length(filelist),") does not match ngrad and is larger then 1\n Please provide each gradient cube in a separate file or one 4D file or use pattern to select.\n")
      if (length(filelist) == ngrad) {
        if (is.null(order)) {
          order <- 1:ngrad
        } else {
          if (length(order) != ngrad)
            stop("readDWIdata: Length of order vector does not match ngrad")
        }
        filelist <- filelist[order]
      }
    } else {
      if (length(filelist) != ngrad)
        stop("readDWIdata: Number of found files does not match ngrad",length(filelist),"\nPlease provide each gradient cube in a separate file.")
      if (is.null(order)) {
        order <- 1:ngrad
      } else {
        if (length(order) != ngrad)
          stop("readDWIdata: Length of order vector does not match ngrad")
      }
      filelist <- filelist[order]
    }
  }
  nfiles <- length(filelist)
  
  ## read all files
  if (verbose) cat("readDWIdata: Start reading data", format(Sys.time()), "\n")
  si <- numeric()
  ddim <- NULL
  first <- TRUE
  i <- 0
  if (!verbose) pb <- txtProgressBar(0, nfiles, style = 3)
  for (ff in filelist) {
    i <- i+1
    if (!verbose) setTxtProgressBar(pb, i)
    if (format == "DICOM") {
      dd <- readDICOMFile(ff, skipSequence = TRUE)
      delta <- c(as.numeric(unlist(strsplit(extractHeader(dd$hdr, "PixelSpacing", FALSE)[1], " "))), extractHeader(dd$hdr, "SliceThickness")[1])
      imageOrientationPatient <- as.numeric(unlist(strsplit(extractHeader(dd$hdr, "ImageOrientationPatient", FALSE)[1], " ")))
      imageOrientationPatient <- matrix(c(imageOrientationPatient, vcrossp(imageOrientationPatient[1:3], imageOrientationPatient[4:6])), 3, 3)
      gradx <- dd$hdr[which((dd$hdr[, 1] == "0019") & (dd$hdr[, 2] == "10BB"))[1], 6]
      grady <- dd$hdr[which((dd$hdr[, 1] == "0019") & (dd$hdr[, 2] == "10BC"))[1], 6]
      gradz <- dd$hdr[which((dd$hdr[, 1] == "0019") & (dd$hdr[, 2] == "10BD"))[1], 6]
      bvalueDCM <- as.numeric(unlist(strsplit(dd$hdr[which((dd$hdr[, 1] == "0043") & (dd$hdr[, 2] == "1039"))[1], 6], " ")))[1]
      if (verbose) cat("diffusion gradient", gradx, grady, gradz, "b-value", bvalueDCM, "\n")
      #      ## WORKAROUND!!
      #      dd$img <- aperm(dd$img, c(2, 1))
      #      ## END WORKAROUND!!
    } else if (format == "NIFTI") {
      dd <- readNIfTI(ff, reorient = FALSE)
      nslice <- dim(dd)[3]
      if (is.null(zind)) zind <- 1:nslice
      delta <- dd@pixdim[2:4]
      imageOrientationPatient <- t(matrix(c(dd@srow_x[1:3]/dd@pixdim[2:4], dd@srow_y[1:3]/dd@pixdim[2:4], dd@srow_z[1:3]/dd@pixdim[2:4]), 3, 3))
    } else if (format == "ANALYZE") {
      dd <- readANALYZE(ff)
      if ( SPM2) dd@.Data <- dd@.Data * dd@funused1
      nslice <- dim(dd)[3]
      if (is.null(zind)) zind <- 1:nslice
      delta <- dd@pixdim[2:4]
      imageOrientationPatient <- diag(3)
      if ( length( dim( dd)) == 4) if ( dim( dd)[ 4] == 1) dim( dd) <-  dim( dd)[ 1:3]
    } else if (format == "AFNI") {
      dd <- readAFNI(ff)
      nslice <- dim(dd)[3]
      if (is.null(zind)) zind <- 1:nslice
      delta <- dd@DELTA
      imageOrientationPatient <- diag(3)
    } 
    ddim <- if (format == "DICOM") c(dim(dd$img)[1:2], nslice, ngrad) else c(dim(dd)[1:2], nslice, ngrad)
    
    if (is.null(voxelext)) {
      if (length(delta) == 3) {
        voxelext <- delta
      } else {
        voxelext <- c(1, 1, 1)
        warning("readDWIdata: Could not find voxel size. Setting default.")
      }
    } else {
      if (length(delta) == 3) {
        if (any(voxelext != delta))
          warning("readDWIdata: Voxel extension", voxelext, "is not match its value in data files:", delta)
      }
    }
    
    if (is.null(rotation)) {
      if (any(imageOrientationPatient != 0)) {
        rotation <- imageOrientationPatient
      } else {
        rotation <- diag(3)
      }
    } else {
      if (any(imageOrientationPatient != rotation)) {
        warning("readDWIdata: Rotation matrices differ: ", imageOrientationPatient)
      }
    }
    
    if (is.null(xind)) xind <- 1:ddim[1]
    if (is.null(yind)) yind <- 1:ddim[2]
    if (format == "DICOM") {
      if (first) { 
        ttt <- dd$img[xind, yind]
        nttt <- dim(ttt)
        si <- numeric(nfiles * prod(nttt))
        dim(si) <- c(nttt, nfiles)
        si[ , , 1] <- ttt
        first <- FALSE
      } else {
        si[ , , i] <- dd$img[xind, yind]
      }
    } else {
      if (length(filelist) > 1) { # list of 3D files
        ## for dti_leipzig data we got length(dim(dd))==4
        if(length(dim(dd))==4&&dim(dd)[4]==1) dim(dd) <- dim(dd)[1:3]
        if (first) { 
          ttt <- dd[xind, yind, zind]
          nttt <- dim(ttt)
          si <- numeric(nfiles * prod(nttt))
          dim(si) <- c(nttt, nfiles)
          si[ , , , 1] <- ttt
          first <- FALSE
        } else {
          si[ , , , i] <- dd[xind, yind, zind]
        }
      } else { # this is a 4D file
        si <- dd[xind, yind, zind,]
      }
    }
  }
  if (!verbose) close(pb)
  dim(si) <- c(length(xind), length(yind), length(zind), ngrad)
  dimsi <- dim(si)
  if (verbose) cat("readDWIdata: Data successfully read", format(Sys.time()), "\n")
  
  # redefine orientation
  xyz <- (orientation)%/%2+1
  swap <- orientation%%2
  if(any(xyz!=1:3)) {
    abc <- 1:3
    abc[xyz] <- abc
    si <- aperm(si,c(abc,4))
    swap[xyz] <- swap
    voxelext[xyz] <- voxelext
    dimsi[xyz] <- dimsi[1:3]
    ddim[xyz] <- ddim[1:3]
    gradient[xyz,] <- gradient
  }
  if(swap[1]==1) {
    si <- si[dimsi[1]:1,,,] 
    gradient[1,] <- -gradient[1,]
  }
  if(swap[2]==1) {
    si <- si[,dimsi[2]:1,,]  
    gradient[2,] <- -gradient[2,]
  }
  if(swap[3]==0) {
    si <- si[,,dimsi[3]:1,]    
    gradient[3,] <- -gradient[3,]
  }
  # orientation set to radiological convention
  
  ## this replaces the content off all voxel with elements <=0 or >maxvalue by 0
  si <- .Fortran("initdata",
                 si = as.double(si),
                 as.integer(dimsi[1]),
                 as.integer(dimsi[2]),
                 as.integer(dimsi[3]),
                 as.integer(dimsi[4]),
                 as.double(maxvalue),
                 PACKAGE = "dti")$si
  if(all(si==as.integer(si))) si <- as.integer(si)
  ##  reduce memory requirements if possible
  dim(si) <- dimsi
  
  ## set level to level*mean  of positive s_0 values
  level <- max(mins0value, level * mean(si[ , , , s0ind][si[ , , , s0ind] > 0]))
  if (verbose) cat("readDWIdata: Create auxiliary statistics",format(Sys.time()), " \n")
  design <- create.designmatrix.dti(gradient)
  
  invisible(new("dtiData",
                call        = args,
                si          = si,
                gradient    = gradient,
                bvalue      = as.vector(bvalue),
                btb         = sweep( design, 2, bvalue, "*"),
                ngrad       = ngrad,
                s0ind       = s0ind,
                replind     = replind(gradient),
                ddim        = as.integer(dim(si)[1:3]),
                ddim0       = as.integer(ddim),
                xind        = xind,
                yind        = yind,
                zind        = zind,
                level       = level,
                sdcoef      = c(1,0,level,quantile(si[,,,s0ind],.98)),
                voxelext    = voxelext,
                orientation = orientation,
                rotation    = rotation,
                source      = paste(dirlist, collapse = "|"))
  )
}

################################################################
#                                                              #
# Section for summary(), print(), show() functions (generic)   #
#                                                              #
################################################################

setMethod("print", "dtiData",
          function(x){
            cat("  Object of class", class(x),"\n")
            cat("  Generated by calls   :\n")
            print(x@call)
            cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
            ns0 <- length(x@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(x@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", x@source, "\n")
            cat("  Slots                :\n")
            print(slotNames(x))
            invisible(NULL)
          })
setMethod("print", "dtiTensor",
          function(x){
            cat("  Object of class", class(x),"\n")
            cat("  Generated by calls   :\n")
            print(x@call)
            cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
            ns0 <- length(x@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(x@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", x@source, "\n")
            cat("  Slots                :\n")
            print(slotNames(x))
            invisible(NULL)
          })

setMethod( "print", "dkiTensor",
           function( x) {
             cat( "  Object of class", class(x), "\n")
             cat( "  Generated by calls   :\n")
             print( x@call)
             cat( "  Dimension            :", paste( x@ddim, collapse = "x"), "\n")
             ns0 <- length( x@s0ind)
             cat( "  Number of S0 images  :", paste( ns0, collapse = "x"), "\n")
             cat( "  Number of Gradients  :", paste( x@ngrad - ns0, collapse = "x"), "\n")
             cat( "  Source-Filename      :", x@source, "\n")
             cat( "  Slots                :\n")
             print( slotNames(x))
             invisible( NULL)
           })

setMethod("print", "dwiMixtensor",
          function(x){
            cat("  Object of class", class(x),"\n")
            cat("  Generated by calls   :\n")
            print(x@call)
            cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
            ns0 <- length(x@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(x@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", x@source, "\n")
            cat("  Naximal number od mixture components:",max(x@order),"\n") 
            cat("  Slots                :\n")
            print(slotNames(x))
            invisible(NULL)
          })
setMethod("print", "dwiQball",
          function(x){
            cat("  Object of class", class(x),"\n")
            cat("  Generated by calls   :\n")
            print(x@call)
            cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
            ns0 <- length(x@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(x@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", x@source, "\n")
            cat("  Kind                 :", paste(x@what, collapse="x"), "\n")
            cat("  Order                :", paste(x@order, collapse="x"), "\n")
            cat("  Slots                :\n")
            print(slotNames(x))
            invisible(NULL)
          })
setMethod("print","dtiIndices",
          function(x){
            cat("  Object of class", class(x),"\n")
            cat("  Generated by calls   :\n")
            print(x@call)
            cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
            ns0 <- length(x@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(x@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", x@source, "\n")
            cat("  Slots                :\n")
            print(slotNames(x))
            invisible(NULL)
          })

setMethod("print","dkiIndices",
          function(x){
            cat("  Object of class", class(x),"\n")
            cat("  Generated by calls   :\n")
            print(x@call)
            cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
            ns0 <- length(x@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(x@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", x@source, "\n")
            cat("  Slots                :\n")
            print(slotNames(x))
            invisible(NULL)
          })

setMethod("print","dwiFiber",
          function(x){
            cat("  Object of class", class(x),"\n")
            cat("  Generated by calls   :\n")
            print(x@call)
            cat("  Dimension            :", paste(x@ddim, collapse="x"), "\n")
            ns0 <- length(x@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(x@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", x@source, "\n")
            cat("  Minimum FA  :", x@minfa, "\n")
            cat("  Maximum angle :", x@maxangle , "\n")
            cat("  Slots                :\n")
            print(slotNames(x))
            invisible(NULL)
          })

setMethod("show", "dtiData",
          function(object){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls   :\n")
            print(object@call)
            cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", object@source, "\n")
            cat("  Slots                :\n")
            print(slotNames(object))
            invisible(NULL)
          })
setMethod("show", "dtiTensor",
          function(object){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls   :\n")
            print(object@call)
            cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", object@source, "\n")
            cat("  Slots                :\n")
            print(slotNames(object))
            invisible(NULL)
          })

setMethod( "show", "dkiTensor",
           function( object) {
             cat( "  Object of class", class( object), "\n")
             cat( "  Generated by calls   :\n")
             print( object@call)
             cat( "  Dimension            :", paste( object@ddim, collapse = "x"), "\n")
             ns0 <- length( object@s0ind)
             cat( "  Number of S0 images  :", paste( ns0, collapse = "x"), "\n")
             cat( "  Number of Gradients  :", paste( object@ngrad - ns0, collapse = "x"), "\n")
             cat( "  Source-Filename      :", object@source, "\n")
             cat( "  Slots                :\n")
             print( slotNames(object))
             invisible( NULL)
           })

setMethod("show", "dwiMixtensor",
          function(object){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls   :\n")
            print(object@call)
            cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", object@source, "\n")
            cat("  Naximal number od mixture components:",max(object@order),"\n") 
            cat("  Slots                :\n")
            print(slotNames(object))
            invisible(NULL)
          })
setMethod("show", "dtiIndices",
          function(object){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls   :\n")
            print(object@call)
            cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", object@source, "\n")
            cat("  Slots                :\n")
            print(slotNames(object))
            invisible(NULL)
          })

setMethod("show", "dkiIndices",
          function(object){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls   :\n")
            print(object@call)
            cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", object@source, "\n")
            cat("  Slots                :\n")
            print(slotNames(object))
            invisible(NULL)
          })

setMethod("show","dwiFiber",
          function(object){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls   :\n")
            print(object@call)
            cat("  Dimension            :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images  :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients  :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Source-Filename      :", object@source, "\n")
            cat("  Minimum FA  :", object@minfa, "\n")
            cat("  Maximum angle :", object@maxangle , "\n")
            cat("  Slots                :\n")
            print(slotNames(object))
            invisible(NULL)
          })

setMethod("summary", "dtiData",
          function(object, ...){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls    :\n")
            print(object@call)
            cat("  Source-Filename       :", object@source, "\n")
            cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images   :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients   :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Voxel extensions      :", paste(signif(object@voxelext,3), collapse="x"), "\n")
            cat("  Index of S0-Images    :", paste(object@s0ind, collapse="x"), "\n")
            cat("  Quantiles of S0-values:","\n")
            print(signif(quantile(object@si[,,,object@s0ind],...),3))
            cat("  Mean S0-value         :", paste(z <- signif(mean(object@si[,,,object@s0ind]),3),collapse="x"), "\n")
            cat("  Threshold for mask    :", paste(signif(object@level,3),collapse="x"), "\n")
            cat("\n")
            invisible(NULL)
          })
setMethod("summary", "dtiTensor",
          function(object, ...){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls    :\n")
            print(object@call)
            cat("  Source-Filename       :", object@source, "\n")
            cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
            ngrad <- object@ngrad
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images   :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients   :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Voxel extensions      :", paste(signif(object@voxelext,3), collapse="x"), "\n")
            cat("  Quantiles of S0-values:","\n")
            print(signif(quantile(object@th0,...),3))
            cat("  Mean S0-value         :", paste(z <- signif(mean(object@th0),3),collapse="x"), "\n")
            cat("  Voxel in mask         :", paste(sum(object@mask), collapse="x"), "\n")
            cat("  Spatial smoothness    :", paste(signif(object@bw,3), collapse="x"), "\n")
            cat("  mean variance         :", paste(signif(mean(object@sigma[object@mask]),3), collapse="x"), "\n")
            penBIC <- log(ngrad-ns0)/(ngrad-ns0)*6
            cat("  BIC         :", paste(signif(mean(log(object@sigma[object@mask]))+penBIC,3), collapse="x"), "\n")
            cat("  hmax                  :", paste(object@hmax, collapse="x"), "\n")
            if(length(object@outlier)>0) cat("  Number of outliers    :", paste(length(object@outlier), collapse="x"), "\n")
            cat("\n")
            invisible(NULL)
          })

setMethod( "summary", "dkiTensor",
           function( object, ...) {
             cat( "  Object of class", class( object), "\n")
             cat( "  Generated by calls    :\n")
             print( object@call)
             cat( "  Source-Filename       :", object@source, "\n")
             cat( "  Dimension             :", paste( object@ddim, collapse = "x"), "\n")
             ngrad <- object@ngrad
             ns0 <- length( object@s0ind)
             cat( "  Number of S0 images   :", paste( ns0, collapse = "x"), "\n")
             cat( "  Number of Gradients   :", paste( object@ngrad - ns0, collapse = "x"), "\n")
             cat( "  Voxel extensions      :", paste( object@voxelext, collapse = "x"), "\n")
             cat( "  Quantiles of S0-values:", "\n")
             print( signif( quantile( object@th0, ...), 3))
             cat( "  Mean S0-value         :", paste( z <- signif( mean( object@th0), 3), collapse = "x"), "\n")
             cat( "  Voxel in mask         :", paste( sum( object@mask), collapse = "x"), "\n")
             cat( "  Spatial smoothness    :", paste( signif( object@bw, 3), collapse = "x"), "\n")
             cat( "  mean variance         :", paste( signif( mean( object@sigma[ object@mask]), 3), collapse = "x"), "\n")
             penBIC <- log( ngrad - ns0) / ( ngrad - ns0) * 6
             cat( "  BIC                   :", paste( signif( mean( log( object@sigma[ object@mask])) + penBIC, 3), collapse = "x"), "\n")
             cat( "  hmax                  :", paste( object@hmax, collapse = "x"), "\n")
             if( length( object@outlier) > 0) cat( "  Number of outliers    :", paste( length( object@outlier), collapse = "x"), "\n")
             cat( "\n")
             invisible( NULL)
           })

setMethod("summary", "dwiMixtensor",
          function(object, ...){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls    :\n")
            print(object@call)
            cat("  Source-Filename       :", object@source, "\n")
            cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
            ngrad <- object@ngrad
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images   :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients   :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Voxel extensions      :", paste(signif(object@voxelext,3), collapse="x"), "\n")
            cat("  Quantiles of S0-values:","\n")
            print(signif(quantile(object@th0,...),3))
            cat("  Mean S0-value         :", paste(z <- signif(mean(object@th0),3),collapse="x"), "\n")
            cat("  Voxel in mask         :", paste(sum(object@mask), collapse="x"), "\n")
            cat("  Spatial smoothness    :", paste(signif(object@bw,3), collapse="x"), "\n")
            cat("  mean variance         :", paste(signif(mean(object@sigma[object@mask]),3), collapse="x"), "\n")
            penBIC <- log(ngrad-ns0)/(ngrad-ns0)*(1+2*object@order[object@mask])
            cat("  BIC         :", paste(signif(mean(log(object@sigma[object@mask])+penBIC),3), collapse="x"), "\n")
            cat("  hmax                  :", paste(object@hmax, collapse="x"), "\n")
            if(length(object@outlier)>0) cat("  Number of outliers    :", paste(length(object@outlier), collapse="x"), "\n")
            nofmc <- table(object@order[object@mask])
            cat("  Numbers of mixture components:") 
            cat(paste(names(nofmc),": ",nofmc,"  ",sep=""))
            cat("\n\n")
            invisible(NULL)
          })
setMethod("summary", "dwiQball",
          function(object, ...){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls    :\n")
            print(object@call)
            cat("  Source-Filename       :", object@source, "\n")
            cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images   :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients   :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Voxel extensions      :", paste(signif(object@voxelext,3), collapse="x"), "\n")
            cat("  Kind                  :", paste(object@what, collapse="x"), "\n")
            ord <- object@order
            cat("  Order                 :", paste(ord, collapse="x"), "\n")
            cat("  Quantiles of S0-values:","\n")
            print(signif(quantile(object@th0,...),3))
            cat("  Mean S0-value         :", paste(z <- signif(mean(object@th0),3),collapse="x"), "\n")
            cat("  Voxel in mask         :", paste(sum(object@mask), collapse="x"), "\n")
            cat("  Spatial smoothness    :", paste(signif(object@bw,3), collapse="x"), "\n")
            cat("  mean variance         :", paste(signif(mean(object@sigma[object@mask]),3), collapse="x"), "\n")
            penBIC <- log(object@ngrad-ns0)*(ord+1)*(ord+2)/2/(object@ngrad-ns0)
            cat("  BIC         :", paste(signif(mean(log(object@sigma[object@mask]))+penBIC,3), collapse="x"), "\n")
            cat("  hmax                  :", paste(object@hmax, collapse="x"), "\n")
            if(length(object@outlier)>0) cat("  Number of outliers    :", paste(length(object@outlier), collapse="x"), "\n")
            cat("\n")
            invisible(NULL)
          })
setMethod("summary", "dtiIndices",
          function(object, ...){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls    :\n")
            print(object@call)
            cat("  Source-Filename       :", object@source, "\n")
            cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images   :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients   :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Voxel extensions      :", paste(signif(object@voxelext,3), collapse="x"), "\n")
            cat("  Percentage of zero values      :",paste(signif(mean(object@fa==0)*100,3), "%",collapse="x"), "\n")
            cat("  Quantiles of positive FA-values:","\n")
            print(signif(quantile(object@fa[object@fa>0],...),3))
            cat("  Quantiles of positive GA-values:","\n")
            print(signif(quantile(object@ga[object@ga>0],...),3))
            cat("  Quantiles of positive MD-values:","\n")
            print(signif(quantile(object@md[object@md>0],...),3))
            cat("\n")
            invisible(NULL)
          })

setMethod("summary", "dkiIndices",
          function(object, ...){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls    :\n")
            print(object@call)
            cat("  Source-Filename       :", object@source, "\n")
            cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images   :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients   :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Voxel extensions      :", paste(signif(object@voxelext,3), collapse="x"), "\n")
            cat("  Percentage of zero values      :",paste(signif(mean(object@fa==0)*100,3), "%",collapse="x"), "\n")
            cat("  Quantiles of positive FA-values:","\n")
            print(signif(quantile(object@fa[object@fa>0],...),3))
            cat("  Quantiles of positive GA-values:","\n")
            print(signif(quantile(object@ga[object@ga>0],...),3))
            cat("  Quantiles of positive MD-values:","\n")
            print(signif(quantile(object@md[object@md>0],...),3))
            cat("\n")
            invisible(NULL)
          })

setMethod("summary","dwiFiber",
          function(object){
            cat("  Object of class", class(object),"\n")
            cat("  Generated by calls    :\n")
            print(object@call)
            cat("  Source-Filename       :", object@source, "\n")
            cat("  Dimension             :", paste(object@ddim, collapse="x"), "\n")
            ns0 <- length(object@s0ind)
            cat("  Number of S0 images   :", paste(ns0, collapse="x"), "\n")
            cat("  Number of Gradients   :", paste(object@ngrad-ns0, collapse="x"), "\n")
            cat("  Voxel extensions      :", paste(signif(object@voxelext,3), collapse="x"), "\n")
            cat("  Minimum FA  :", object@minfa, "\n")
            cat("  Maximum angle :", object@maxangle , "\n")
            cat("  Number of fibers :", length(object@startind), "\n")
            cat("  Quantiles of fiber lengths:\n")
            print(quantile(diff(c(object@startind,dim(object@fibers)[1]+1)-2)))
            cat("  Total number of line segments :", dim(object@fibers)[1]-length(object@startind),"\n")
            #  linesegments in one fiber = length(fiber - 1)
            cat("\n")
            invisible(NULL)
          })

################################################################
#                                                              #
# Section for interface functions                              #
# like: tensor2medinria()                                      #
#                                                              #
################################################################

tensor2medinria <- function(obj, filename, xind=NULL, yind=NULL, zind=NULL) {
  
  if (is.null(xind)) xind <- 1:obj@ddim[1]
  if (is.null(yind)) yind <- 1:obj@ddim[2]
  if (is.null(zind)) zind <- 1:obj@ddim[3]
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  
  D <- aperm( obj@D, c( 2:4, 1))[ xind, yind, zind, c( 1, 2, 4, 3, 5, 6)]
  dim(D) <- c( length(xind), length(yind), length(zind), 1, 6)
  nim <- nifti(D,
               dim_ = c( 5, length(xind), length(yind), length(zind), 1, 6, 1, 1),
               pixdim = c( -1, obj@voxelext[1:3], 1, 1, 0, 0),
               intent_code = 1007,
               datatype = 16,
               bitpix = 32, ## must correspond to datatype
               sclslope = 1,
               xyztunits = "\002", # ???
               qform = 1,
               sform = 1,
               quatern_d = 1,
               srow_x = c( -obj@voxelext[1], 0, 0, 0),
               srow_y = c( 0, obj@voxelext[2], 0, 0),
               srow_z = c( 0, 0, obj@voxelext[3], 0)
  )
  
  writeNIfTI( nim, filename)
}

medinria2tensor <- function(filename) {
  args <- sys.call() 
  data <- readNIfTI(filename, reorient = FALSE)
  
  invisible(new("dtiTensor",
                call  = list(args),
                D     = aperm(data, c( 5, 1:4))[ c( 1, 2, 4, 3, 5, 6), , , , , drop = TRUE],
                sigma = array(0, dim(data)[1:3]),
                scorr = array(0, c( 1, 1, 1)),
                bw    = rep( 0, 3),
                mask  = array(TRUE, dim(data)[1:3]),
                method = "unknown",
                hmax  = 1,
                th0   = array(0, dim = dim(data)[1:3]),
                gradient = matrix(0,1,1),
                btb   = matrix(0,1,1),
                ngrad = as.integer(0), # = dim(btb)[2]
                s0ind = as.integer(0),
                ddim  = dim(data)[1:3],
                ddim0 = dim(data)[1:3],
                xind  = 1:dim(data)[1],
                yind  = 1:dim(data)[2],
                zind  = 1:dim(data)[3],
                voxelext = data@pixdim[2:4],
                orientation = as.integer(c(0,2,5)),
                rotation = diag(3),
                scale = 1,
                source= "unknown")
  )
  
}
