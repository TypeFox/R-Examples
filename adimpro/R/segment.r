segment <- function(object, level=0.5, delta=0, thresh= 3,
                    fov=NULL,channel=0, hmax=4, aws=TRUE, varmodel=NULL,
                    ladjust=1.25 , xind = NULL,
                    yind = NULL, wghts=c(0.299, 0.587, 0.114, 0),
                    scorr=TRUE, lkern="Triangle",
                    plateau=NULL, homogen=TRUE, earlystop=TRUE, demo=FALSE, select=FALSE, sext=1.4, connected=FALSE,
                    graph=FALSE, max.pixel=4.e2,compress=TRUE) {
  #
  #
  #          Auxilary functions
  #
  IQRdiff <- function(data) IQR(diff(data))/1.908
  #
  #   sequence of factors for lambda obtained by new propagation condition 
  #
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
#
#  first extract
#
  if(level < 0 || level >65535) return(warning("Illegal value of level"))
  if(level-delta < 0 || level+delta >65535) return(warning("Illegal value of delta"))
  if(level + delta <= 1) {
     level <- 65535 * level
     delta <- 65535 * delta
  }
  img <- extract.image(object)
  dimg <- dimg0 <- dim(img)
  if(length(dimg)>2) {
     if(channel > dimg[3] || channel==0){
        dim(img) <- c(dimg[1]*dimg[2],dimg[3])
        img <- wghts[1:dimg[3]] %*% t(img)
        dim(img) <- object$dim
  } else  {
        img <- img[,,channel]
  }
  }
  imgtype <- object$type 
  if(!is.null(object$call)) { 
    warning("argument object is result of awsimage or awspimage. You should know what you do.")
  }
  args <- match.call()
#  determine parameters for segmentation
  if(select){
     graph <- TRUE
cat("select region inside the interesting structure\n")
     timg <- extract.image(clip.image(make.image(img)))
     level <- median(timg)
     delta  <- delta+sext*IQR(timg)/1.34898
cat("Specified level: ",level,"  delta: ",delta,"\n")
#  thats ext1 times estimated standard deviation 
  } 
  if(is.null(varmodel)) {
    varmodel <- if(object$gamma) "Constant" else "Linear"
  }
  bcf <- c(-.4724,1.1969,.22167,.48691,.13731,-1.0648)
# coefficients for  bias correction for spatial correlation
  estvar <- toupper(varmodel) %in% c("CONSTANT","LINEAR","QUADRATIC")
  hpre <- 2.
  hvest <- 3.5
  if(estvar) {
    dlw<-(2*trunc(hpre)+1)
    nvarpar <- switch(toupper(varmodel),CONSTANT=1,LINEAR=2,QUADRATIC=3,1)
  } else {
    nvarpar <- 1
  }
  #
  #   Check image type
  #
  spcorr <- numeric(2)
  h0 <- c(0,0)
  #
    if(is.null(xind)) xind <- 1:dimg[1]
    if(is.null(yind)) yind <- 1:dimg[2]
    n1 <- length(xind)
    n2 <- length(yind)
    n <- n1*n2
  dimg <- c(n1,n2)
  #
  #     set approriate defaults
  #
#  ladjust <- max(1,ladjust)
  lkern <- switch(lkern,
                  Triangle=1,
                  Quadratic=2,
                  Cubic=3,
                  Plateau=4,
                  1)
  if (is.null(hmax)) hmax <- 4
  wghts <- wghts/sum(wghts)
  if (aws) lambda <- ladjust*switch(imgtype,greyscale=16,rgb=7) else lambda <- 1e50
  #
  #      in case of colored noise get the corresponding bandwidth (for Gaussian kernel)
  #
  sigma2 <- max(.01,IQRdiff(img)^2)
  cat("Estimated variance (assuming independence): ", signif(sigma2/65635^2,4),"\n")
  if (!estvar) {
    if (length(sigma2)==1) {
      #   homoskedastic Gaussian case
    } 
  }
    # set the support of the statistical kernel to (0,1), set spmin 
    spmin <- plateau
    if(is.null(spmin)) spmin <- .25
# determine maximum volume (variance reduction)
  maxvol <- getvofh2(hmax,lkern)
  if(is.null(fov)) fov <- maxvol
  kstar <- as.integer(log(maxvol)/log(1.25))  
  if(aws){ 
     k <- if(estvar) 6 else 1 
     }
  else {
    cat("No adaptation method specified. Calculate kernel estimate with bandwidth hmax.\n")
    k <- kstar
  }
  if (demo && !graph) graph <- TRUE
  if(graph){
    oldpar <- par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
    on.exit(par(oldpar))
    if(exists(".adimpro")&&!is.null(.adimpro$xsize)) max.x <- trunc(.adimpro$xsize/3.3)  else max.x <- 500
    if(exists(".adimpro")&&!is.null(.adimpro$xsize)) max.y <- trunc(.adimpro$xsize/1.2)  else max.y <- 1000
#  specify settings for show.image depending on geometry (mfrow) and maximal screen size
  }
  # now check which procedure is appropriate
  #
  #    Initialize  list for theta
  #
  bi <- rep(1,n)
  theta <- img[xind,yind]
  segment <- matrix(0,n1,n2)
  bi0 <- 1
  #
  #  if varmodel specified prepare for initial variance estimation
  #
    coef <- numeric(nvarpar)
    coef[1] <- sigma2
    varest <- matrix(sigma2,n1,n2)
    vobj <- list(coef=coef,meanvar=sigma2)
    imgq995 <- quantile(img[xind,yind],.995)
  if(scorr){
    twohp1 <- 2*trunc(hpre)+1
    pretheta <- .Fortran("awsimg",
                         as.integer(img[xind,yind]),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(1),
                         as.double(hpre),
                         theta=integer(prod(dimg)),
                         bi=double(n1*n2),
                         as.integer(lkern),
                         double(twohp1*twohp1),# array for location weights
                         double(1),DUP=FALSE,
                         PACKAGE="adimpro")$theta
    dim(pretheta) <- dimg
    spchcorr <- .Fortran("estcorr",
                         as.double(img[xind,yind]-pretheta),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(1),
                         scorr=double(2),
                         double(1),
                         as.double(hpre),
                         DUP=FALSE,
                         PACKAGE="adimpro")[c("scorr")]
    spcorr <- spchcorr$scorr
    srh <- sqrt(hpre) 
    spcorr <- pmin(.9,spcorr+
                          bcf[1]/srh+bcf[2]/hpre+
                          bcf[3]*spcorr/srh+bcf[4]*spcorr/hpre+
                          bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hpre)
    #  bias correction for spatial correlation
    cat("Estimated spatial correlation in channel ",signif(spcorr,2),"\n")
  } else {
    spcorr <- rep(0,2)
  }
  ###
  ###              gridded   2D
  ###
  lambda0 <- 1e50
  #
  #   run single steps to display intermediate results
  #
   total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
#
#
#   Start main loop
#
#
   while (k<=kstar) {
    hakt0 <- geth2(1,10,lkern,1.25^(k-1),1e-4)
    hakt <- geth2(1,10,lkern,1.25^k,1e-4)
    twohp1 <- 2*trunc(hakt)+1
    spcorr <- pmax(spcorr,0)
    if(any(spcorr>0)) {
      h0<-numeric(length(spcorr))
      for(i in 1:length(h0))
        h0[i]<-geth.gauss(spcorr[i])
      if(length(h0)<2) h0<-rep(h0[1],2)
      # cat("Corresponding bandwiths for specified correlation:",h0,"\n")
    }
    if (any(spcorr>0)) {
      lcorr <- Spatialvar.gauss(hakt0,h0)/
        Spatialvar.gauss(h0,1e-5)/Spatialvar.gauss(hakt0,1e-5)
      # Correction C(h0,hakt) for spatial correlation depends on h^{(k-1)} 
      lambda0 <-lambda0*lcorr
      if(varmodel=="None") 
        lambda0 <- lambda0*Varcor.gauss(h0)
    } 
      zobj <- .Fortran("segment",
                       as.integer(img[xind,yind]),
                       as.double(level),
                       as.double(delta),
                       as.integer(n1),
                       as.integer(n2),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.integer(theta),
                       as.double(vobj$coef),
                       as.integer(nvarpar),
                       as.double(vobj$meanvar),
                       bi=as.double(bi),
                       double(n1*n2),# array for inverse variances
                       theta=integer(n1*n2),
                       as.integer(lkern),
                       as.double(spmin),
                       double(twohp1*twohp1),# array for location weights
                       pvalues=double(n1*n2),# array for pvalues
                       segment=as.integer(segment),# array for segment (-1,0,1)
                       as.double(thresh),
                       as.double(fov),
                       varest=as.double(varest),
                       DUP=FALSE,
                       PACKAGE="adimpro")[c("bi","theta","hakt","pvalues","segment","varest")]
    varest <- pmin(.01,zobj$varest)
    theta <- zobj$theta
    bi <- zobj$bi
    pvalues <- zobj$pvalues
    segment <- zobj$segment
    dim(theta) <- dim(pvalues) <- dim(segment) <- dim(bi) <- c(n1,n2)
    rm(zobj)
    gc()
    if (graph) {
      graphobj <- make.image(img[xind,yind],compress=FALSE)
      show.image(graphobj,max.x=max.x,max.y=max.y,xaxt="n",yaxt="n")
      title("Observed Image")
      graphobj <- make.image(theta,compress=FALSE)
      show.image(graphobj,max.x=max.x,max.y=max.y,xaxt="n",yaxt="n")
      title(paste("Reconstruction  h=",signif(hakt,3)))
      graphobj <- make.image(32767.5*(1+segment),compress=FALSE)
      show.image(graphobj,max.x=max.x,max.y=max.y,xaxt="n",yaxt="n")
      title(paste("Segmentation  h=",signif(hakt,3)))
      graphobj$img <- matrix(as.integer(65534*bi/max(bi)),n1,n2)
      show.image(graphobj,max.x=max.x,max.y=max.y,xaxt="n",yaxt="n")
      title(paste("Adaptation (rel. weights):",signif(mean(bi)/max(bi),3)))
      rm(graphobj)
      gc()
    }
    cat("Bandwidth",signif(hakt,3)," Progress",signif(total[k],2)*100,"% ")
    cat("segmented:",sum(segment==-1),sum(segment==0),sum(segment==1))
    cat("\n")
    cat("range of bi",range(bi),"\n")
    if (scorr) {
      #
      #   Estimate Correlations  (keep old estimates until hmax > hpre)
      #
      if(hakt > hpre & hakt < hvest){
      spchcorr <- .Fortran("estcorr",
                           as.double(img[xind,yind] - theta),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(1),
                           scorr=double(2),
                           double(1),
                           DUP=FALSE,
                           PACKAGE="adimpro")["scorr"]
    spcorr <- spchcorr$scorr
    srh <- sqrt(hakt) 
    spcorr <- pmin(.9,spcorr+
                         bcf[1]/srh+bcf[2]/hakt+
                         bcf[3]*spcorr/srh+bcf[4]*spcorr/hakt+
                         bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hakt)
      #  bias correction for spatial correlation
      cat("Estimated spatial correlation:",signif(spcorr,2),"\n")
      } else {
            spcorr <- pmin(.9,spchcorr$scorr+
                         bcf[1]/srh+bcf[2]/hpre+
                         bcf[3]*spchcorr$scorr/srh+bcf[4]*spchcorr$scorr/hpre+
                         bcf[5]*spchcorr$scorr^2/srh+bcf[6]*spchcorr$scorr^2/hpre)
      #  bias correction for spatial correlation
      }
    } else {
      spcorr <- rep(0,2)
    }
    if (estvar & hakt < hvest) {
      #
      #   Create new variance estimate
      #
      vobj <- .Fortran(switch(toupper(varmodel),CONSTANT="esigmac",LINEAR="esigmal",QUADRATIC="esigmaq"),
                       as.integer(img[xind,yind]),
                       as.integer(n1*n2),
                       as.integer(1),
                       as.integer(theta),
                       as.double(bi),
                       as.integer(imgq995),
                       coef=double(nvarpar),
                       meanvar=double(1),
                       DUP=FALSE,
                       PACKAGE="adimpro")[c("coef","meanvar")]
      cat("Estimated mean variance",signif(vobj$meanvar/65635^2,3),"\n")
    }
    if (demo) readline("Press return")
    lambda0 <- lambda
    k <- k+1
  }
#
#   END main loop
#
#
#   check for connectivity if specified
#
if(connected) {
   show.image(make.image(segment), main = "Identify center  of segmented region by left mouse click")
cat("Identify center  of segmented region by left mouse click\n")
        coord <- locator(1, type = "p")
    segment <- matrix(.Fortran("connect1",
                               segment=as.integer(segment),
                               as.integer(n1),
                               as.integer(n2),
                               as.integer(coord$x),
                               as.integer(coord$y),
                               integer(n1*n2),
                               integer(n1*n2),
                               logical(n1*n2),
                               DUP=FALSE,
                               PACKAGE="adimpro")$segment-1,n1,n2)
}
#
#
#
  ###                                                                       
  ###            end cases                                                  
  ###                                 .....................................
  if(graph) par(oldpar)
  object <- make.image(segment)
  object$xind <- xind
  object$yind <- yind
  object$hsegm <- hakt
  object$level <- level
  object$delta <- delta
  object$thresh <- thresh
  object$call <- args
  invisible(if(compress) compress.image(object) else object)
}








