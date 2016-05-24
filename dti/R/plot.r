################################################################
#                                                              #
# Section for plot() functions                                 #
#                                                              #
################################################################

setMethod("plot", "dwi", function(x, y, ...) cat("No implementation for class dwi\n"))

##############

setMethod("plot", "dtiData", function(x, y,slice=1, gradient=NULL, view= "axial", show=TRUE, density=FALSE, xind=NULL, yind=NULL, zind=NULL, mar=c(3,3,3,.3), mgp=c(2,1,0), ...) {
  if(is.null(x@si)) cat("No dwi data yet")
  maxsi <- max(x@si)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if(is.null(gradient)) gradient <- x@s0ind[1]
  if(gradient<1||gradient>x@ngrad) {
    warning("gradient number out of range, show s0 image")
    gradient <- x@s0ind[1]
  }
  if(density) { 
    z <- density(x@si[xind,yind,zind,gradient])
    if(show) {
      plot(z,main="Density of S0-values") 
      lines(c(x@level,x@level),c(0,max(z$y)),col=2)
    }
    return(invisible(z))
  }
  #adimpro <- require(adimpro)
  if (view == "sagittal") {
    if(slice<1||slice>x@ddim[1]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[1]%/%2
    }
    img <- x@si[slice,yind,zind,gradient]
  } else if (view == "coronal") {
    if(slice<1||slice>x@ddim[2]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[2]%/%2
    }
    img <- x@si[xind,slice,zind,gradient]
  } else {
    if(slice<1||slice>x@ddim[3]) {
      warning("slice number out of range, show central slice")
      slice <- x@ddim[3]%/%2
    }
    img <- x@si[xind,yind,slice,gradient]
  }
  oldpar <- par(mar=mar,mgp=mgp, ...)
  #if(adimpro) {
  img <- make.image(65535*img/maxsi)
  if(show) show.image(img,...)
  #} else if(show) {
  #  image(img,...)
  #}
  par(oldpar)
  invisible(img)
})

##############

setMethod("plot", "dtiTensor", function(x, y, slice=1, view="axial", quant=0, minfa=NULL, contrast.enh=1,what="fa", qrange=c(.01,.99),xind=NULL,yind=NULL,zind=NULL, mar=c(2,2,2,.2),mgp=c(2,1,0),...) {
  if(is.null(x@D)) cat("No diffusion tensor yet")
  #adimpro <- require(adimpro)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if (view == "sagittal") {
    D <- x@D[,slice,yind,zind]
    mask <- x@mask[slice,yind,zind]
    n1 <- length(yind)
    n2 <- length(zind)
  } else if (view == "coronal") {
    D <- x@D[,xind,slice,zind]
    mask <- x@mask[xind,slice,zind]
    n1 <- length(xind)
    n2 <- length(zind)
  } else {
    D <- x@D[,xind,yind,slice]
    mask <- x@mask[xind,yind,slice]
    n1 <- length(xind)
    n2 <- length(yind)
  }
  if(what=="GA"){
    z <- .Fortran("dti2Dga",
                  as.double(D),
                  as.integer(n1*n2),
                  as.logical(mask),
                  fa=double(n1*n2),
                  md=double(n1*n2),
                  andir=double(3*n1*n2),
                  PACKAGE="dti")[c("fa","md","andir")]
  } else {
    z <- .Fortran("dti2Dfa",
                  as.double(D),
                  as.integer(n1*n2),
                  as.logical(mask),
                  fa=double(n1*n2),
                  md=double(n1*n2),
                  andir=double(3*n1*n2),
                  PACKAGE="dti")[c("fa","md","andir")]
  }
  oldpar <- par(mfrow=c(3,3),mar=mar,mgp=mgp,...)
  #  now draw information to graphical device
  on.exit(par(oldpar))
  img<-D[1,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  show.image(make.image(65535*img/max(img)))
  title(paste("Dxx: mean",signif(mean(D[mask]),3),"max",signif(max(D[1,,][mask]),3)))
  img<-D[2,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  img <- (img-rg[1])/(rg[2]-rg[1])
  show.image(make.image(img))
  title(paste("Dxy: min",signif(min(D[2,,][mask]),3),"max",signif(max(D[2,,][mask]),3)))
  img<-D[3,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  img <- (img-rg[1])/(rg[2]-rg[1])
  show.image(make.image(img))
  title(paste("Dxz: min",signif(min(D[3,,][mask]),3),"max",signif(max(D[3,,][mask]),3)))
  show.image(make.image(matrix(z$fa,n1,n2)))
  if(what=="GA"){
    title(paste("Anisotropy Index  (GA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
  } else {
    title(paste("Fractional Anisotropy (FA)  range:",signif(min(z$fa[mask]),3),"-",
                signif(max(z$fa[mask]),3)))
  }
  img<-D[4,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  show.image(make.image(65535*img/max(img)))
  title(paste("Dyy: min",signif(min(D[4,,][mask]),3),"max",signif(max(D[4,,][mask]),3)))
  img<-D[5,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  img <- (img-rg[1])/(rg[2]-rg[1])
  show.image(make.image(img))
  title(paste("Dyz: min",signif(min(D[5,,][mask]),3),"max",signif(max(D[5,,][mask]),3)))
  andir.image(matrix(z$fa,n1,n2),array(z$andir,c(3,n1,n2)),quant=quant,minfa=minfa)
  title(paste("Anisotropy directions"))
  img <- matrix(z$md,n1,n2)
  show.image(make.image(65535*img/max(img)))
  if(what=="GA"){
    title(paste("Mean log diffusivity   range:",signif(min(z$md[mask]),3),"-",
                signif(max(z$md[mask]),3)))
  } else {
    title(paste("Mean diffusivity   range:",signif(min(z$md[mask]),3),"-",
                signif(max(z$md[mask]),3)))
  }
  img<-D[6,,]
  rg<-quantile(img,qrange)
  img[img>rg[2]]<-rg[2]
  img[img<rg[1]]<-rg[1]
  show.image(make.image(65535*img/max(img)))
  title(paste("Dzz: min",signif(min(D[6,,][mask]),3),"max",signif(max(D[6,,][mask]),3)))
  invisible(NULL)
})

##############

setMethod("plot", "dwiMixtensor", function(x, y, slice=1, view="axial", what="fa", minfa=NULL, identify=FALSE,  xind=NULL,yind=NULL,zind=NULL, mar=c(2,2,2,.2),mgp=c(2,1,0),...) {
  #adimpro <- require(adimpro)
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if (view == "sagittal") {
    x <- x[slice,yind,zind]
  } else if (view == "coronal") {
    x <- x[xind,slice,zind]
  } else {
    x <- x[xind,yind,slice]
  }
  what <- tolower(what)
  stats <- extract(x,what)
  oldpar <- par(mfrow=c(1,length(what)),mar=mar,mgp=mgp,...)
  on.exit(par(oldpar))
  if("fa" %in% what){
    fa <- drop(stats$fa)
    if(!is.null(minfa)) fa[fa<minfa] <- 0
    show.image(img <- make.image(65535*fa))
    title(paste("Slice",slice,"FA"))
  }
  if("order" %in% what){
    order <- drop(stats$order)
    show.image(img <- make.image(65535*order/max(order)))
    title(paste("Slice",slice,"Order of mixture (Maximum=",max(order),")"))
  }
  if("eorder" %in% what){
    eorder <- drop(stats$eorder)
    show.image(img <- make.image(65535*eorder/max(eorder)))
    title(paste("Slice",slice,"Effective order of mixture (Maximum=",signif(max(eorder),2),")"))
  }
  if("ev" %in% what){
    ev <- drop(stats$ev[1,,,])
    show.image(img <- make.image(65535*ev/max(ev)))
    title(paste("Slice",slice,"Maximal Eigenvalue (Maximum=",signif(max(ev),3),")"))
  }
  if(identify){
    xind<-(1:x@ddim[1])
    yind<-(1:x@ddim[2])
    zind<-(1:x@ddim[3])
    img <- extract.image(img)
    image(1:dim(img)[1],1:dim(img)[2],img,col=grey((0:255)/255))
    identifyFA(view,slice,xind,yind,zind)
  } else {
    par(oldpar)
    invisible(img)
  }
})

##############

setMethod("plot", "dtiIndices", function(x, y, slice=1, view= "axial", method=1, quant=0, minfa=NULL, show=TRUE, identify=FALSE, density=FALSE, contrast.enh=1,what="fa",xind=NULL,yind=NULL,zind=NULL, mar=c(3,3,3,.3),mgp=c(2,1,0), ...) {
  what <- tolower(what)
  if(is.null(x@fa)) cat("No anisotropy index yet")
  if(!(method %in% 1:6)) {
    warning("method out of range, reset to 1")
    method <- 1
  }
  if(is.null(xind)) xind<-(1:x@ddim[1])
  if(is.null(yind)) yind<-(1:x@ddim[2])
  if(is.null(zind)) zind<-(1:x@ddim[3])
  if(density) { 
    x <- x[xind,yind,zind]
    z <- density(if(what=="fa") x@fa[x@fa>0] else x@ga[x@ga>0])
    if(show) {
      plot(z,main=paste("Density of positive",what,"-values")) 
    }
    return(invisible(z))
  }
  #adimpro <- require(adimpro)
  oldpar <- par(mar=mar,mgp=mgp, ...)
  if (view == "sagittal") {
    anindex <- if(what=="ga") tanh(x@ga[slice,yind,zind]) else x@fa[slice,yind,zind]
    if (method == 3) {
      andirection <- x@bary[,slice,yind,zind]
    } else {
      andirection <- x@andir[,slice,yind,zind]
    }
  } else if (view == "coronal") {
    anindex <- if(what=="ga") tanh(x@ga[xind,slice,zind]) else x@fa[xind,slice,zind]
    if (method == 3) {
      andirection <- x@bary[,xind,slice,zind]
    } else {
      andirection <- x@andir[,xind,slice,zind]
    }
  } else {
    anindex <- if(what=="ga") tanh(x@ga[xind,yind,slice]) else x@fa[xind,yind,slice]
    if (method == 3) {
      andirection <- x@bary[,xind,yind,slice]
    } else {
      andirection <- x@andir[,xind,yind,slice]
    }
  }
  anindex[anindex>1] <- 0
  anindex[anindex<0] <- 0
  if ((method==1) || (method==2) || (method==4)) {
    if (contrast.enh < 1 && contrast.enh > 0) anindex <- pmin(anindex/contrast.enh, 1)
    if(is.null(minfa)) minfa <- quantile(anindex,quant,na.rm=TRUE)
    if (diff(range(anindex,na.rm=TRUE)) == 0) minfa <- 0
    if(method==1) {
      andirection[1,,] <- abs(andirection[1,,])
      andirection[2,,] <- abs(andirection[2,,])
      andirection[3,,] <- abs(andirection[3,,])
    } else if (method==2) {
      ind<-andirection[3,,]<0
      dim(andirection) <- c(3,prod(dim(ind)))
      andirection[,ind] <- - andirection[,ind]
      andirection[1,] <- (1+andirection[1,])/2
      andirection[2,] <- (1+andirection[2,])/2
      andirection <- sweep(andirection,2,sqrt(apply(andirection^2,2,sum)),"/")
      dim(andirection) <- c(3,dim(ind))
    } else {
      andirection[1,,] <- andirection[1,,]^2
      andirection[2,,] <- andirection[2,,]^2
      andirection[3,,] <- andirection[3,,]^2
    }
    andirection <- aperm(andirection,c(2,3,1))
    andirection <- andirection*as.vector(anindex)*as.numeric(anindex>minfa)
    #    if(adimpro) {
    andirection[is.na(andirection)] <- 0
    andirection <- make.image(andirection,gammatype="ITU")
    if(show) show.image(andirection,...)
    if(identify){
      identifyFA(view,slice,xind,yind,zind)
    } else {
      par(oldpar)
      invisible(andirection)
    }
    #    } else if(show) {
    #      dim(anindex) <- dim(andirection)[2:3]
    #      image(anindex,...)
    #      if(identify){
    #         identifyFA(view,slice,xind,yind,zind)
    #      } else {
    #        par(oldpar)
    #        invisible(NULL)
    #      }
    #    }
  } else if (method==3) {
    #    if(adimpro) {
    andirection[is.na(andirection)] <- 0
    bary <- make.image(aperm(andirection,c(2,3,1)))
    if(show) show.image(bary,...)
    if(identify){
      identifyFA(view,slice,xind,yind,zind)
    } else {
      par(oldpar)
      invisible(bary)
    }
    #    } else if(show) {
    #      image(andirection[1,,],...)
    #      if(identify){
    #         identifyFA(view,slice,xind,yind,zind)
    #      } else {
    #         par(oldpar)
    #         invisible(NULL)
    #      }
    #    }
  } else if (method==5) {
    #    if(adimpro) {
    andirection[is.na(andirection)] <- 0
    img.hsi.data <- array(0,dim=c(dim(andirection)[2:3],3))
    img.hsi.data[,,1] <- atan2(andirection[2,,],andirection[1,,])
    img.hsi.data[,,1] <- img.hsi.data[,,1] + pi*(img.hsi.data[,,1]<0)
    img.hsi.data[,,2] <- abs(acos(andirection[3,,]))
    img.hsi.data[,,3] <- anindex
    img.hsi <- make.image(img.hsi.data,gammatype="ITU",xmode="HSI")
    if(show) show.image(img.hsi,...)
    if(identify){
      identifyFA(view,slice,xind,yind,zind)
    } else {
      par(oldpar)
      invisible(img.hsi)
    }
    #    } else if(show) {
    #      image(andirection[1,,],...)
    #      if(identify){
    #         identifyFA(view,slice,xind,yind,zind)
    #      } else {
    #         par(oldpar)
    #         invisible(NULL)
    #      }
    #    }
  } else if (method==6) {
    data("colqFA", envir = environment())
    #    if (adimpro) {
    colqFA <- col2rgb(colqFA)/255
    img.data <- array(0, dim=c(dim(anindex), 3))
    for (i in 1:dim(anindex)[1]) { # i dont like for loops in R!
      for (j in 1:dim(anindex)[2]) {
        img.data[i,j,] <- colqFA[, 255 * anindex[i, j] + 1] 
      }
    }
    img <- make.image(img.data, gammatype="ITU")
    if(show) show.image(img, ...)
    if(identify){
      identifyFA(view,slice,xind,yind,zind)
    } else {
      par(oldpar)
      invisible(img)
    }
    #    } else if (show) {
    #      image(anindex, col=colqFA)
    #      if(identify){
    #         identifyFA(view,slice,xind,yind,zind)
    #      } else {
    #         par(oldpar)
    #         invisible(NULL)
    #      }
    #    }
  }
})

setMethod("plot", "dwiFiber", function(x, y, ...) {
  plot(density(diff(c(x@startind,dim(x@fibers)[1]+1))-1), ...)
  title("Density of fiber lenghts")
  invisible(NULL)
})

setMethod("plot", "dkiIndices", function(x, 
                                         y, 
                                         slice = 1, 
                                         #view = "axial", 
                                         what = c("md", "fa", "mk", "mk2"),
                                         #method = 1, 
                                         #quant = 0, 
                                         #minfa = NULL, 
                                         #show = TRUE, 
                                         #identify = FALSE, 
                                         #density = FALSE, 
                                         #contrast.enh = 1,
                                         #xind = NULL,
                                         #yind = NULL,
                                         #zind = NULL, 
                                         #mar = c(3, 3, 3, .3),
                                         #mgp = c(2, 1, 0), 
                                         ...) {
  
  switch(what,
         md = image(x@md[, , slice], col=grey(0:255/255), ...),
         fa = image(x@fa[, , slice], col=grey(0:255/255), ...),
         mk = image(x@mk[, , slice], col=grey(0:255/255), ...),
         mk2 = image(x@mk2[, , slice], col=grey(0:255/255), ...))
  
  
  
})
