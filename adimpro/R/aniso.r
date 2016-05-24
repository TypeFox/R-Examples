imganiso2D <- function(x,satexp=.25,g=3,rho=0){
# x should be an array of dimension (3,n1,n2)
require(adimpro)
if(class(x)=="adimpro") x <- tensor2D(x,g=g,rho=rho)
x[c(1,3),,] <- x[c(1,3),,]*(1+1.e-8)+rho
p <- x[1,,]+x[3,,]
q <- x[1,,]*x[3,,]-x[2,,]^2
z <- p^2/4-q
z[z<0] <- 0
# just to avoid numerical -0's
ew <- p/2 + sqrt(z)
ew2 <- p/2 - sqrt(z)
theta <- acos(x[2,,]/sqrt((ew-x[1,,])^2+x[2,,]^2))/pi
theta[is.na(theta)] <- 0
img <- array(0,c(dim(x)[2:3],3))
img[,,1] <- theta
img[,,2] <- (ew/max(ew))^satexp
img[,,3] <- log(ew/ew2)
invisible(make.image(img,xmode="hsi"))
}



tensor2D <- function(x,g=0,rho=0){
if(class(x)=="adimpro") {
   x <- extract.image(x)
}
storage.mode(x) <- "numeric"
ddim <- dim(x)
n1 <- ddim[1]
n2 <- ddim[2]
if(length(ddim)==2){
dx <- rbind(x,x[n1,])-rbind(x[1,],x)
dy <- cbind(x,x[,n2])-cbind(x[,1],x)
dxx <- (dx[-1,]^2+dx[-n1-1,]^2)/2
dyy <- (dy[,-1]^2+dy[,-n2-1]^2)/2
dxy <- (dx[-1,]*dy[,-1]+dx[-n1-1,]*dy[,-1]+dx[-1,]*dy[,-n2-1]+dx[-n1,]*dy[,-n2-1])/4
z<-aperm(array(c(dxx,dxy,dyy),c(ddim,3)),c(3,1,2))
} else {
z <- array(0,c(n1,n2,3))
for(i in 1:3) {
dx <- rbind(x[,,i],x[n1,,i])-rbind(x[1,,i],x[,,i])
dy <- cbind(x[,,i],x[,n2,i])-cbind(x[,1,i],x[,,i])
dxx <- (dx[-1,]^2+dx[-n1-1,]^2)/2
dyy <- (dy[,-1]^2+dy[,-n2-1]^2)/2
dxy <- (dx[-1,]*dy[,-1]+dx[-n1-1,]*dy[,-1]+dx[-1,]*dy[,-n2-1]+dx[-n1,]*dy[,-n2-1])/4
z <- z+c(dxx,dxy,dyy)
}
z <- aperm(z/3,c(3,1,2))
}
if(g>0){
z <- array(.Fortran("sm2Dtens",
                     as.double(z),
                     as.integer(n1),
                     as.integer(n2),
                     as.double(g),
                     as.double(rho),
                     zhat=double(3*n1*n2),
                     PACKAGE="adimpro")$zhat,dim(z))
z[is.na(z)] <- 0
}
z
}

awsaniso <- function (object, hmax=4, g=3, rho=0, aws=TRUE, varmodel=NULL,
                      ladjust=1.0,  xind = NULL, 
                      yind = NULL, wghts=c(1,1,1,1), scorr=TRUE, lkern="Triangle", 
                      demo=FALSE, graph=FALSE, satexp=0.25, max.pixel=4.e2,clip=FALSE,compress=TRUE) {

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
  object <- switch(object$type,
                   "hsi" = hsi2rgb(object),
                   "yuv" = yuv2rgb(object),
                   "yiq" = yiq2rgb(object),
                   "xyz" = xyz2rgb(object),
                   object)
  if(object$type=="RAW") stop("RAW-image, will be implemented in function awsraw later ")
  if(!is.null(object$call)) { 
    warning("argument object is result of awsimage or awspimage. You should know what you do.")
  }
  args <- match.call()
  if(clip) {
        object <- clip.image(object,xind,yind,compress=FALSE)
        xind <- NULL
        yind <- NULL
        }
  if(object$compressed) object <- decompress.image(object)
  if(is.null(varmodel)) {
    varmodel <- if(object$gamma) "Constant" else "Linear"
  }
  bcf <- c(-.4724,1.1969,.22167,.48691,.13731,-1.0648)
# coefficients for  bias correction for spatial correlation
  estvar <- toupper(varmodel) %in% c("CONSTANT","LINEAR")
  hpre <- 2.
  if(estvar) {
    dlw<-(2*trunc(hpre)+1)
    nvarpar <- switch(varmodel,Constant=1,Linear=2,1)
  }
  #
  #   Check image type
  #
  dimg <- dimg0 <- dim(object$img)
  imgtype <- object$type 
  dv <- switch(imgtype,rgb=dimg[3],greyscale=1)
  if(length(wghts)<dv) wghts <- c(wghts,rep(1,dv-length(wghts))) else wghts <- pmin(.1,wghts) 
  wghts <- switch(imgtype, greyscale=1, rgb = wghts[1:dv])
  spcorr <- matrix(0,2,dv)
  h0 <- c(0,0)
    if(is.null(xind)) xind <- 1:dimg[1]
    if(is.null(yind)) yind <- 1:dimg[2]
    n1 <- length(xind)
    n2 <- length(yind)
    n <- n1*n2
  dimg <- switch(imgtype,greyscale=c(n1,n2),rgb=c(n1,n2,dv))
  #
  #     set approriate defaults
  #
  ladjust <- max(1,ladjust)
  lkern <- switch(lkern,
                  Triangle=1,
                  Quadratic=2,
                  Cubic=3,
                  Plateau=4,
                  1)
  if (is.null(hmax)) hmax <- 4
  wghts <- wghts/sum(wghts)
  dgf <- sum(wghts)^2/sum(wghts^2)
  if (aws) lambda <- ladjust*switch(imgtype, "greyscale" = 8.1, "rgb" = 4.4, 4.4) else lambda <- 1e50
  #
  #      in case of colored noise get the corresponding bandwidth (for Gaussian kernel)
  #
  sigma2 <- pmax(0.01,switch(imgtype,
                   "greyscale" = IQRdiff(object$img)^2,
                   "rgb" = apply(object$img,3,IQRdiff)^2,
                   IQRdiff(object$img)^2))
  cat("Estimated variance (assuming independence): ", signif(sigma2/65635^2,4),"\n")
  if (!estvar) {
    if (length(sigma2)==1) {
      #   homoskedastic Gaussian case
      if(dv>1) sigma2 <- rep(sigma2,1)
      lambda <- lambda*sigma2 
    } else if (length(sigma2)==dv) {
      wghts <- wghts[1:dv]/sigma2
    } 
  }
    # set the support of the statistical kernel to (0,1), set spmin and spmax
    lambda <- 2*lambda
    spmin <- 0
    spmax <- 1
# determine maximum volume (variance reduction)
  maxvol <- getvofh2(hmax,lkern)
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
    graphobj0 <- object[-(1:length(object))[names(object)=="img"]]
    graphobj0$dim <- c(n1,n2)
    if(exists(".adimpro")&&!is.null(.adimpro$xsize)) max.x <- trunc(.adimpro$xsize/3.3)  else max.x <- 500
    if(exists(".adimpro")&&!is.null(.adimpro$xsize)) max.y <- trunc(.adimpro$xsize/1.2)  else max.y <- 1000
#  specify settings for show.image depending on geometry (mfrow) and maximal screen size
  }
  # now check which procedure is appropriate
  #
  #    Initialize  list for theta
  #
  bi <- rep(1,n)
  theta <- switch(imgtype,greyscale=object$img[xind,yind],rgb=object$img[xind,yind,])
  #
  #  if varmodel specified prepare for initial variance estimation
  #
  if(estvar){
    coef <- matrix(0,nvarpar,dv)
    coef[1,] <- sigma2
    vobj <- list(coef=coef,meanvar=sigma2)
    imgq995 <- switch(imgtype,greyscale=quantile(object$img[xind,yind],.995),
                      rgb=apply(object$img[xind,yind,],3,quantile,.995))
  }
  if(scorr){
    twohp1 <- 2*trunc(hpre)+1
    pretheta <- .Fortran("awsimg",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         as.double(hpre),
                         theta=integer(prod(dimg)),
                         double(n1*n2),
                         as.integer(lkern),
                         double(twohp1*twohp1),# array for location weights
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")$theta
    dim(pretheta) <- dimg
    spchcorr <- .Fortran("estcorr",
                         as.double(switch(imgtype,
                                          greyscale=object$img[xind,yind],
                                          rgb=object$img[xind,yind,])-pretheta),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         scorr=double(2*dv),
                         chcorr=double(max(1,dv*(dv-1)/2)),
                         as.double(hpre),
                         DUP=FALSE,
                         PACKAGE="adimpro")[c("scorr","chcorr")]
    spcorr <- spchcorr$scorr
    srh <- sqrt(hpre) 
    spcorr <- matrix(pmin(.9,spcorr+
                          bcf[1]/srh+bcf[2]/hpre+
                          bcf[3]*spcorr/srh+bcf[4]*spcorr/hpre+
                          bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hpre),2,dv)
    #  bias correction for spatial correlation
    chcorr <- spchcorr$chcorr
    for(i in 1:dv) cat("Estimated spatial correlation in channel ",i,":",signif(spcorr[,i],2),"\n")
    if(dv>1) cat("Estimated correlation between channels (1,2):",signif(chcorr[1],2),
        " (1,3):",signif(chcorr[2],2),
        " (2,3):",signif(chcorr[3],2),"\n")
  } else {
    spcorr <- matrix(0,2,dv)
    chcorr <- numeric(max(1,dv*(dv-1)/2))
  }
  #
  #         fix values of the image in inactiv pixel
  #
  ###
  ###              gridded   2D
  ###
  lambda0 <- 1e50
  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  #
  #   run single steps to display intermediate results
  #
   while (k<=kstar) {
    hakt0 <- geth2(1,10,lkern,1.25^(k-1),1e-4)
    hakt <- geth2(1,10,lkern,1.25^k,1e-4)
    twohp1 <- 2*trunc(hakt)+1
    spcorr <- pmax(apply(spcorr,1,mean),0)
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
    if(dv==1) dim(theta) <- c(n1,n2) else dim(theta) <- c(n1,n2,dv)
    anisoimg <- tensor2D(theta,g=g,rho=max(rho,1e-6)*mean(vobj$meanvar)/bi)
    anisoimg[1,,] <- anisoimg[1,,]+max(rho,1e-6)*mean(vobj$meanvar)/bi
    anisoimg[3,,] <- anisoimg[3,,]+max(rho,1e-6)*mean(vobj$meanvar)/bi
    hakt0 <- hakt
      if(estvar){
        varcase <- 1
        zobj <- .Fortran("aniawsv",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         as.double(anisoimg),
                         as.double(vobj$coef),
                         as.integer(nvarpar),
                         as.double(vobj$meanvar),
                         as.double(chcorr),
                         hakt=as.double(hakt),
                         as.double(lambda0),
                         as.integer(theta),
                         bi=as.double(bi),
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         as.integer(1),
                         as.double(spmin),		       
                         as.double(spmax),
                         as.double(sqrt(wghts)),
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")[c("bi","theta","hakt")]
      } else {
        # all other cases
        varcase <- 3
        zobj <- .Fortran("aniawsim",
                         as.integer(switch(imgtype,
                                           greyscale=object$img[xind,yind],
                                           rgb=object$img[xind,yind,])),
                         as.integer(n1),
                         as.integer(n2),
                         as.integer(dv),
                         as.double(anisoimg),
                         hakt=as.double(hakt),
                         as.double(lambda0),
                         as.integer(theta),
                         bi=as.double(bi),
                         theta=integer(prod(dimg)),
                         as.integer(lkern),
                         as.integer(1),
                         as.double(spmin),
                         as.double(spmax),
                         as.double(wghts),
                         double(dv),DUP=FALSE,
                         PACKAGE="adimpro")[c("bi","theta","hakt")]
      } 
    theta <- zobj$theta
    bi <- zobj$bi
    rm(zobj)
    gc()
    dim(bi) <- dimg[1:2]
    if (graph) {
      graphobj <- graphobj0
      class(graphobj) <- "adimpro"
      graphobj$img <- switch(imgtype,
                             greyscale=object$img[xind,yind],
                             rgb=object$img[xind,yind,])
      show.image(graphobj,max.x=max.x,max.y=max.y,xaxt="n",yaxt="n")
      title("Observed Image")
      graphobj$img <- array(as.integer(theta),dimg)
      show.image(graphobj,max.x=max.x,max.y=max.y,xaxt="n",yaxt="n")
      title(paste("Reconstruction  h=",signif(hakt,3)))
      show.image(imganiso2D(anisoimg,satexp=satexp),max.x=max.x,max.y=max.y,xaxt="n",yaxt="n")
      title(paste("Estimated anisotropy  g=",signif(g,3)))
      graphobj$img <- matrix(as.integer(65534*bi/max(bi)),n1,n2)
      graphobj$type <- "greyscale"
      graphobj$gamma <- FALSE
      show.image(graphobj,max.x=max.x,max.y=max.y,xaxt="n",yaxt="n")
      title(paste("Adaptation (rel. weights):",signif(mean(bi)/max(bi),3)))
      rm(graphobj)
      gc()
    }
    cat("Bandwidth",signif(hakt,3)," Progress",signif(total[k],2)*100,"% \n")
    if (scorr) {
      #
      #   Estimate Correlations  (keep old estimates until hmax > hpre)
      #
      if(hakt > hpre){
      spchcorr <- .Fortran("estcorr",
                           as.double(switch(imgtype,
                                            greyscale=object$img[xind,yind],
                                            rgb=object$img[xind,yind,]) - theta),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(dv),
                           scorr=double(2*dv),
                           chcorr=double(max(1,dv*(dv-1)/2)),
                           DUP=FALSE,
                           PACKAGE="adimpro")[c("scorr","chcorr")]
    spcorr <- spchcorr$scorr
#    spcorr <- matrix(pmin(.9,0.8817*spcorr+0.231/hakt+6.018*spcorr/hakt^2+
#                            1.753*spcorr^2/hakt-10.622*spcorr^2/hakt^2),2,dv)
    srh <- sqrt(hakt) 
    spcorr <- matrix(pmin(.9,spcorr+
                         bcf[1]/srh+bcf[2]/hakt+
                         bcf[3]*spcorr/srh+bcf[4]*spcorr/hakt+
                         bcf[5]*spcorr^2/srh+bcf[6]*spcorr^2/hakt),2,dv)
      #  bias correction for spatial correlation
      chcorr <- spchcorr$chcorr
      for(i in 1:dv) cat("Estimated spatial correlation in channel ",i,":",signif(spcorr[,i],2),"\n")
      if(dv>1) cat("Estimated correlation between channels (1,2):",signif(chcorr[1],2),
                   " (1,3):",signif(chcorr[2],2),
                   " (2,3):",signif(chcorr[3],2),"\n")
      } else {
           spcorr <-   matrix(pmin(.9,spchcorr$scorr+
                         bcf[1]/srh+bcf[2]/hpre+
                         bcf[3]*spchcorr$scorr/srh+bcf[4]*spchcorr$scorr/hpre+
                         bcf[5]*spchcorr$scorr^2/srh+bcf[6]*spchcorr$scorr^2/hpre),2,dv)         
# spcorr <- matrix(pmin(.9,.06+0.9*spchcorr$scorr+1.108/hpre),2,dv)
      #  bias correction for spatial correlation
         chcorr <- spchcorr$chcorr
      }
    } else {
      spcorr <- matrix(0,2,dv)
      chcorr <- numeric(max(1,dv*(dv-1)/2))
    }
    if (estvar) {
      #
      #   Create new variance estimate
      #
      vobj <- .Fortran(switch(varmodel,Constant="esigmac",Linear="esigmal"),
                       as.integer(switch(imgtype,
                                         greyscale=object$img[xind,yind],
                                         rgb=object$img[xind,yind,])),
                       as.integer(n1*n2),
                       as.integer(dv),
                       as.integer(theta),
                       as.double(bi),
                       as.integer(imgq995),
                       coef=double(nvarpar*dv),
                       meanvar=double(dv),
                       DUP=FALSE,
                       PACKAGE="adimpro")[c("coef","meanvar")]
      dim(vobj$coef) <- c(nvarpar,dv)
       for(i in 1:dv){
      if(any(spcorr[,i]>.1) & vobj$meanvar[i]<sigma2[i]) {
            vobj$coef[,i] <- vobj$coef[,i]*sigma2[i]/vobj$meanvar[i]
            vobj$meanvar[i] <- sigma2[i]
      }
     }
      cat("Estimated mean variance",signif(vobj$meanvar/65635^2,3),"\n")
    }
    k <- k + 1
    if (demo) readline("Press return")
    lambda0 <- lambda*ladjust
  }
  ###                                                                       
  ###            end cases                                                  
  ###                                 .....................................
  if(graph) par(oldpar)
  if(imgtype=="rgb") {
    object$img[xind,yind,] <- theta
  } else if(imgtype=="greyscale") {
    object$img[xind,yind] <- theta
  }
  #  if(dv==1) dim(img) <- dim(img)[1:2]
  ni <- array(1,dimg0[1:2])
  ni[xind,yind] <- bi
  object$ni <- ni
  object$ni0 <- max(ni)
  object$hmax <- hakt
  object$call <- args
  if(estvar) {
    if(varmodel=="Linear") {
       vobj$coef[1,] <- vobj$coef[1,]/65535^2
       vobj$coef[2,] <- vobj$coef[2,]/65535
    } else vobj$coef <- vobj$coef/65535^2
    object$varcoef <- vobj$coef
    object$wghts <- vobj$wghts
  }
  if(scorr) {
    object$scorr <- spcorr
    object$chcorr <- chcorr
  }
  invisible(if(compress) compress.image(object) else object)
}
