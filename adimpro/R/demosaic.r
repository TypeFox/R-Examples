awsraw <- function (object, hmax=4, aws=TRUE, wb=c(1,1,1), cspace="Adobe", ladjust=1.,
                    maxrange=TRUE,lkern="Triangle", 
                    graph=FALSE, max.pixel=4.e2, compress=TRUE) {
#
#  demosaicing using median of 3x3 s_{ij}'s in the smoothing algorithm instead of
#  a separate demosaicing algorithm
#
  #
  #          Auxilary functions
  #
  IQRdiff <- function(data) IQR(diff(data))/1.908

  #####################################################################################
  ###    
  ###    function body
  ###    
  ###    first check arguments and initialize
  ###
  #####################################################################################    
  if(!check.adimpro(object)) {
    cat(" Consistency check for argument object failed (see warnings). object is returned.\"n")
    return(invisible(object)) 
  }
  if(object$type!="RAW") stop("object does not contain RAW sensor data, 
                    please read the image by read.raw(filename,type=''RAW'')")
  if(!is.numeric(wb)) wb <- c(1,1,1)
  if(length(wb)<3) wb <- rep(wb[1],3)
# 
#  detect type of Bayer mask
# 
  args <- match.call()
  if(object$compressed) object <- decompress.image(object)
  varmodel <- "Linear"
  nvarpar <- 2
  #
  #   Check image type
  #
  dimg <- dimg0 <- dim(object$img)
  dv <- 3
  spcorr <- matrix(0,2,dv)
  h0 <- c(0,0)
  n1 <- dimg[1]
  n2 <- dimg[2]
  n <- n1*n2
  #
  #     set approriate defaults
  #
  lkern <- switch(lkern,
                  Triangle=1,
                  Quadratic=2,
                  Cubic=3,
                  Plateau=4,
                  1)
#  if (is.null(hmax)) hmax <- 4
  if(hmax<sqrt(2)) hmax<-sqrt(2)
  dgf <- 3
  if (aws) lambda <- 15*ladjust else lambda <- 1e50
#  provides a max. rel. decrease in weights of 9% for hmax <= 10, ladjust=1
#  compared to nonadaptive smoothing for a homogeneus RAW-image
  cat("using lambda=",lambda,"\n")
  maxvol <- getvofh2(hmax/sqrt(2),lkern)
  kstar <- as.integer(log(maxvol)/log(1.25))  
  k <- 1
  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  #  set the support of the statistical kernel to (0,1), set spmin
    spmin <- .25
  if(graph){
    oldpar <- par(mfrow=c(1,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
    on.exit(par(oldpar))
    graphobj <- object[-(1:length(object))[names(object)%in%c("img","type")]]
    class(graphobj) <- "adimpro"
  }
  # now check which procedure is appropriate
  #
  #    Initialize  list for theta
  #
  bi0 <- 1
  #
  #   extract bayer mask
  #
  bayer <- switch(extract.info(object),RGGB=1,GRBG=2,BGGR=3,GBRG=4)
  if(extract.info(object,"Isize")!=extract.info(object,"Osize")) bayer <- bayer+1
  bayer <- (bayer-1)%%4+1
  if(!is.null(object$rotate)) {
      bayer <- object$rotate+bayer 
      bayer <- (bayer-1)%%4+1
  }
#  if(any(wb!=1)){
# White balance if specified
#     object$img <- matrix(.Fortran("wbalance",
#                              sensor=as.integer(object$img),
#                              as.integer(n1),
#                              as.integer(n2),
#                              as.double(wb),
#                              as.integer(bayer),
#                              DUP=FALSE,
#                              PACKAGE="adimpro")$sensor,n1,n2)
#  }
  if(maxrange){
     minimg <- min(object$img)
     rangeimg <- max(object$img)-minimg
     object$img <- matrix(as.integer((object$img-minimg)/rangeimg*65535),n1,n2)
  }
  dimg <- dimg0 <- dim(object$img)
  n1 <- dimg[1]
  n2 <- dimg[2]
  n <- n1*n2
  out.cam <- cam2rgbmat(object, cspace=cspace)
  #
  #   prepare for initial variance estimation
  #
  coef <- matrix(0,nvarpar,dv)
  vobj <- list(coef=coef,meanvar=c(0,0,0))
  imghom <- .Fortran("dhomogen",
                    as.integer(object$img),
                    as.integer(object$dim[1]),
                    as.integer(object$dim[2]),
                    imghom=integer(n1*n2),
                    as.integer(bayer),
                    DUP=FALSE,
                    PACKAGE="adimpro")$imghom
  medianhom <- median(imghom)
  indnothom <- imghom > 3*medianhom
  rm(imghom)
  gc()
  cat("Using ", n1*n2-sum(indnothom)," of ",n1*n2," pixel for variance estimates\n")
#
#   this should keep 90% of the observations under homogeneity
#                     
    sensorhat0 <- .Fortran("smsens0",
                    as.integer(object$img),
                    shat=integer(n1*n2),
                    bi=double(n1*n2),
                    as.integer(object$dim[1]),
                    as.integer(object$dim[2]),
                    as.integer(bayer),
                    DUP=FALSE,
                    PACKAGE="adimpro")[c("shat","bi")]

  vobj <- .Fortran("senvar",
                    as.integer(object$img),
                    as.integer(object$dim[1]),
                    as.integer(object$dim[2]),
                    as.integer(sensorhat0$shat),
                    as.double(sensorhat0$bi),
                    as.integer(bayer),
                    coef=double(6),
                    meanvar=double(3),
                    as.logical(indnothom),
                    DUP=FALSE,
                    PACKAGE="adimpro")[c("coef","meanvar")]
      dim(vobj$coef) <- c(2,3)
      cat("Estimated mean variance",signif(vobj$meanvar/65635^2,3),"\n")
      cat("Estimated constant variance term",signif(vobj$coef[1,]/65635^2,3),"\n")
      cat("Estimated linear  variance  term",signif(vobj$coef[2,]/65635,3),"\n")
  gc()
  #
  #         fix values of the image in inactiv pixel
  #
  ###
  ###              gridded   2D
  ###
  #
  #   run single steps to display intermediate results
  #
  zobj <- list(bi=matrix(1,n1,n2),shat=object$img)
  cimg <- array(.Fortran("demmed4",
                         as.integer(object$img),
                         cimg=integer((n1-2)*(n2-2)*3),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n1-2),
                         as.integer(n2-2),
                         as.integer(bayer),
                         DUP=FALSE,
                         PACKAGE="adimpro")$cimg,c((n1-2),(n2-2),3))
   while (k<=kstar) {
    hakt0 <- geth2(1,10,lkern,1.25^(k-1),1e-4)*sqrt(2)
    hakt <- geth2(1,10,lkern,1.25^k,1e-4)*sqrt(2)
    twohp1 <- 2*trunc(hakt)+1
    twohp3 <- twohp1+2
    zobj <- .Fortran("smsensor",
                         as.integer(object$img),
                         shat=integer(n1*n2),
                         as.integer(cimg),# thats the color image
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n1-2),
                         as.integer(n2-2),
                         as.integer(bayer),
                         as.double(vobj$coef),
                         as.double(vobj$meanvar),
                         hakt=as.double(hakt),
                         as.double(lambda),
                         bi=as.double(zobj$bi),
                         as.integer(lkern),
                         as.double(spmin), 
                         double(twohp1*twohp1),# array for location weights
                         DUP=FALSE,
                         PACKAGE="adimpro")[c("bi","shat")]
  dim(zobj$bi)  <- c(n1,n2)
  dim(zobj$shat) <- c(n1,n2)
  bi0 <- max(zobj$bi)
  cimg <- array(.Fortran("demmed4b",
                         as.integer(zobj$shat),
                         cimg=as.integer(cimg),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(n1-2),
                         as.integer(n2-2),
                         as.integer(bayer),
                         DUP=FALSE,
                         PACKAGE="adimpro")$cimg,c((n1-2),(n2-2),3))
    gc()
    if(any(wb!=1)){
      cimg[,,1] <- as.integer(pmin(65535,wb[1]*cimg[,,1]))
      cimg[,,2] <- as.integer(pmin(65535,wb[2]*cimg[,,2]))
      cimg[,,3] <- as.integer(pmin(65535,wb[3]*cimg[,,3]))
    }
    if (graph) {
      graphobj$type <- "rgb"
      graphobj$dim <- c(n1-2,n2-2)
      graphobj$img <- array(.Fortran("cam2rgb",
                                      as.integer(cimg),
                                      as.integer((n1-2)*(n2-2)),
                                      as.double(out.cam),
                                      theta=integer((n1-2)*(n2-2)*3),
                                      DUP=FALSE,
                                      PACKAGE="adimpro")$theta,c(dimg-2,3))
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title(paste("Reconstruction  h=",signif(hakt,3)))
      graphobj$dim <- c(n1/2,n2/2)
      bi <- zobj$bi[seq(1,n1,2),seq(1,n2,2)]+zobj$bi[seq(2,n1,2),seq(1,n2,2)]+
            zobj$bi[seq(1,n1,2),seq(2,n2,2)]+zobj$bi[seq(2,n1,2),seq(2,n2,2)]
      graphobj$img <- matrix(as.integer(65534*bi/max(bi)),n1/2,n2/2)
      graphobj$type <- "greyscale"
      graphobj$gamma <- FALSE
      show.image(graphobj,max.x=max.pixel,xaxt="n",yaxt="n")
      title(paste("Adaptation (rel. weights):",signif(mean(bi)/max(bi),3)))
      gc()
    }
      #
      #   Create new variance estimate
      #
      if(hakt>3){
      vobj <- .Fortran("senvar",
                    as.integer(object$img),
                    as.integer(object$dim[1]),
                    as.integer(object$dim[2]),
                    as.integer(zobj$shat),
                    as.double(zobj$bi),
                    as.integer(bayer),
                    coef=double(6),
                    meanvar=double(3),
                    as.logical(indnothom),
                    DUP=FALSE,
                    PACKAGE="adimpro")[c("coef","meanvar")]
      dim(vobj$coef) <- c(2,3)
      cat("Estimated mean variance",signif(vobj$meanvar/65635^2,3),"\n")
      cat("Estimated constant variance term",signif(vobj$coef[1,]/65635^2,3),"\n")
      cat("Estimated linear  variance  term",signif(vobj$coef[2,]/65635,3),"\n")
    }
    cat("Bandwidth",signif(hakt,3)," Progress",signif(total[k],2)*100,"% \n")
    k <- k+1
  }
  ###                                                                       
  ###            end cases                                                  
  ###                                 .....................................
  if(graph) par(oldpar)
  object$img <- array(.Fortran("cam2rgb",
                                as.integer(cimg),
                                as.integer((n1-2)*(n2-2)),
                                as.double(out.cam),
                                theta=integer((n1-2)*(n2-2)*3),
                                DUP=FALSE,
                                PACKAGE="adimpro")$theta,c(dimg-2,3))
  object$type <- "rgb"
  object$cspace <- cspace
  object$ni <- zobj$bi/65535
  object$ni0 <- bi0/65535
  object$dim <- c(n1-2,n2-2)
  object$hmax <- hakt
  object$call <- args
       vobj$coef[1,] <- vobj$coef[1,]/65535^2
       vobj$coef[2,] <- vobj$coef[2,]/65535
  object$vcoef <- vobj$coef
  invisible(if(compress) compress.image(object) else object)
}





