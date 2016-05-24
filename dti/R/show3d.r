################################################################
#                                                              #
# Section for show3d() functions (public)                      #
#                                                              #
################################################################

show3d <- function(obj,  ...) cat("3D Visualization not implemented for this class:",class(obj),"\n")

setGeneric("show3d", function(obj,  ...) standardGeneric("show3d"))

setMethod("show3d","dtiData", function( obj, xind=NULL, yind=NULL, zind=NULL,
                                        quant=.8, scale=.4, bgcolor="black", add=FALSE, maxobjects=729, 
                                        what=c("adc","data"), minalpha=1, nn=1, normalize=FALSE, 
                                        box=FALSE, title=FALSE,...){
  ## check what
  what <- tolower(what)
  what <- match.arg(what)
  #if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  cube <- selectCube(xind,yind,zind,obj@ddim,maxobjects)
  xind <- cube$xind
  yind <- cube$yind
  zind <- cube$zind
  n <- cube$n
  n1 <- cube$n1
  n2 <- cube$n2
  n3 <- cube$n3
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  obj <- obj[xind,yind,zind,drop=FALSE]
  vext <- obj@voxelext
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  dim(tmean) <- c(3,n)
  radii <- extract(obj,"sb")$sb
  s0 <- extract(obj,"s0")$s0
  if(length(dim(s0))==4) s0 <- apply(s0,1:3,mean)
  radii <- sweep(radii,1:3,s0,"/")
  if(what=="adc") radii <- array(pmax(0,-log(radii)),dim(radii))
  # avoid using negative ADC's caused by random effects 
  ngrad <- dim(radii)[length(dim(radii))]
  dim(radii) <- c(length(radii)/ngrad,ngrad)
  radii <- t(radii)
  sscale <- scale
  #  if(what=="colorcoded") sscale <- 1
  if(normalize){
    minradii <- apply(radii,2,min)
    maxradii <- apply(radii,2,max)
    radii <- sweep(radii,2,minradii,"-")
    radii <- sweep(radii,2,maxradii-minradii,"/")*sscale
  } else {
    radii <- radii/quantile(apply(radii,2,max),quant)*sscale
  }
  gradient <- obj@gradient[,-obj@s0ind]
  if(!add) {
    open3d()
    par3d(...)
    rgl.bg(color=bgcolor)
  }
  show3dData(radii,gradient,centers=tmean,minalpha=minalpha,...)
  if(box) bbox3d()
  if(is.character(title)) {
    title3d(title,color="white",cex=1.5)
  } else {
    if(title) title3d(switch(what,"data"="observed DWI data","adc"="observed ADC"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(what,"data"="observed diffusion weighted data","adc"="apparent diffusion coefficients from data"),"\n",
      if(normalize) "normalized","\n")
  invisible(rgl.cur())
})
##############

setMethod("show3d","dtiTensor", function(obj, 
                                         xind=NULL, yind=NULL, zind=NULL, method=1, minfa=.3, mask=NULL, fibers=FALSE, 
                                         maxangle=30, level=0, quant=.8, scale=.4, bgcolor="black", add=FALSE, 
                                         subdivide=2, maxobjects=729, what=c("tensor","adc","odf"), odfscale=1, 
                                         minalpha=.25, normalize=NULL, box=FALSE, title=FALSE, ...){
  #if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  ## check what
  what <- tolower(what)
  what <- match.arg(what)
  if(!exists("icosa0")) data("polyeders", envir = environment())

  ## WORKAROUND to make "polyeders" not global variable (for R CMD check) 
  icosa0 <- icosa0
  icosa1 <- icosa1
  icosa2 <- icosa2
  icosa3 <- icosa3
  icosa4 <- icosa4
  ## END WORKAROUND

  if(subdivide<0||subdivide>4) subdivide <- 3
  cube <- selectCube(xind,yind,zind,obj@ddim,maxobjects)
  xind <- cube$xind
  yind <- cube$yind
  zind <- cube$zind
  n <- cube$n
  n1 <- cube$n1
  n2 <- cube$n2
  n3 <- cube$n3
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  mask <- if(is.null(mask)) obj@mask else obj@mask&mask
  mask <- as.vector(mask[xind,yind,zind])
  obj <- obj[xind,yind,zind]
  vext <- obj@voxelext
  n <- prod(obj@ddim) 
  D <- obj@D
  D <- D/max(D)
  dim(D) <- c(6,n)
  mask <- mask & (D[1,]*D[4,]*D[6,]>0)
  tmean <- array(0,c(3,obj@ddim))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  dim(tmean) <- c(3,n)
  z <- extract(obj,what=c("andir","fa"))
  if(minfa>0) mask <- mask&(z$fa>=minfa)
  maxev <- extract(obj,what="evalues",mc.cores=1)$evalues[3,,,,drop=FALSE][mask, drop=FALSE]
  dim(mask) <- NULL
  andir <- matrix(z$andir,3,n)[,mask,drop=FALSE]
  fa <- z$fa[mask, drop=FALSE]
  if(method==1) {
    andir <- abs(andir)
  } else {
    ind<-andir[1,]<0
    andir[,ind] <- - andir[,ind]
    andir[2,] <- (1+andir[2,])/2
    andir[3,] <- (1+andir[3,])/2
  }
  colorvalues <- rgb(andir[1,],andir[2,],andir[3,])
  D <- D[,mask, drop=FALSE]
  tmean <- tmean[,mask, drop=FALSE]
  n <- sum(mask)
  if(is.null(normalize)) normalize <- switch(what,"tensor"=FALSE,"adc"=TRUE,"odf"=FALSE)
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  radii <- .Fortran(switch(what,tensor="ellradii",adc="adcradii",odf="odfradii"),
                    as.double(polyeder$vertices),
                    as.integer(polyeder$nv),
                    as.double(D),
                    as.integer(n),
                    radii=double(n*polyeder$nv),
                    PACKAGE="dti")$radii
  dim(radii) <- c(polyeder$nv,n)
  if(what=="odf") normalize <- FALSE
  if(normalize){
    minradii <- apply(radii,2,min)
    maxradii <- apply(radii,2,max)
    radii <- sweep(radii,2,minradii,"-")
    radii <- sweep(radii,2,maxradii-minradii,"/")*scale
  } else {
    if (what=="odf"){
      #
      #   use a sphere of radius level as baseline for the ODF
      #
      radii <- radii+level
      #
      #   to display results in a form that the volumes are comparable,
      #   need volume elements around vertices that have volume proportional to radii,
      #   i.e. a radial extension of radii^(1/3)
      #   odfscale = 1 corresponds to using radii directly for ODF-values
      #   values inbetween are possible
      #
      radii <- radii^(1/odfscale)
      radii <- radii/quantile(apply(radii,2,max),quant)*scale
    } else {
      radii <- (radii+level)/(quantile(apply(radii,2,max),quant)+level)*scale
    }
  }
  if(!add) {
    open3d()
    par3d(...)
    rgl.bg(color=bgcolor)
  }
  if(what=="odf"){
    show3dODF(radii,polyeder,centers=tmean,minalpha=minalpha,...)
  } else {
    show3dTens(radii,polyeder,centers=tmean,colors=colorvalues,alpha=minalpha+(1-minalpha)*fa)
  }
  if(fibers){
    tracks <- tracking(obj,mask=mask,minfa=minfa,maxangle=maxangle)
    dd <- tracks@fibers
    startind <- tracks@startind
    dd <- expandFibers(dd,startind)$fibers
    rgl.lines(dd[,1]+vext[1]/2,dd[,2]+vext[2]/2,dd[,3]+vext[3]/2,
              color=rgb(abs(dd[,4]),abs(dd[,5]),abs(dd[,6])),
              size=3)
  }  
  if(box) bbox3d()
  if(is.character(title)) {
    title3d(title,color="white",cex=1.5)
  } else {
    if(title) title3d(switch(what,"tensor"="estimated tensors","adc"="estimated ADC (tensor)"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(what,"tensor"="estimated tensors","adc"="apparent diffusion coefficients from estimated tensors"),"\n",
      if(obj@hmax>1) paste("smoothed with hmax=",obj@hmax),if(normalize) "normalized","\n")
  invisible(rgl.cur())
})
setMethod("show3d","dwiMixtensor", function(obj, 
                                            xind=NULL, yind=NULL, zind=NULL, minfa=.3, minorder=1, mineo=1, fibers=FALSE, 
                                            maxangle=30, level=0, quant=.8, scale=.4, bgcolor="black", add=FALSE, 
                                            subdivide=3, maxobjects=729, what=c("odf","axis","both"), odfscale=1, 
                                            minalpha=1, lwd=3, box=FALSE, title=FALSE, ...){
  #if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  ## check what
  what <- tolower(what)
  what <- match.arg(what)
  if(!exists("icosa0")) data("polyeders", envir = environment())

  ## WORKAROUND to make "polyeders" not global variable (for R CMD check) 
  icosa0 <- icosa0
  icosa1 <- icosa1
  icosa2 <- icosa2
  icosa3 <- icosa3
  icosa4 <- icosa4
  ## END WORKAROUND

  if(subdivide<0||subdivide>4) subdivide <- 3
  cube <- selectCube(xind,yind,zind,obj@ddim,maxobjects)
  xind <- cube$xind
  yind <- cube$yind
  zind <- cube$zind
  n <- cube$n
  n1 <- cube$n1
  n2 <- cube$n2
  n3 <- cube$n3
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  obj <- obj[xind,yind,zind]
  mask <- obj@mask
  vext <- obj@voxelext
  scale <- scale*min(vext)
  order <- obj@order
  ev <- obj@ev
  mix <- obj@mix
  maxorder <- dim(mix)[1]
  orient <- obj@orient
  n <- prod(obj@ddim) 
  tmean <- array(0,c(3,obj@ddim))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  if(minfa > 0){
    fa <- extract(obj,"fa")$fa
    mask <- (fa>=minfa)&mask
  }
  if(mineo > 1) mask <- (extract(obj,"eorder")$eorder>=mineo)&mask
  if(minorder > 1) mask <- (order>=minorder)&mask
  dim(ev) <- c(2,n)
  dim(tmean) <- c(3,n)
  dim(orient) <- c(dim(orient)[1:2],n)
  dim(mix) <- c(dim(mix)[1],n)
  if(minfa>0||mineo>1||minorder>1){
    indpos <- (1:n)[mask]
    ev <- ev[,indpos]
    tmean <- tmean[,indpos]
    orient <- orient[,,indpos]
    mix <- mix[,indpos]
    order <- order[indpos]
    fa <- fa[indpos]
    n <- length(indpos)
  }
  gc()
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  if(what %in% c("odf","both")){
    radii <- .Fortran("mixtradi",
                      as.double(polyeder$vertices),
                      as.integer(polyeder$nv),
                      as.double(ev),
                      as.double(orient),
                      as.double(mix),
                      as.integer(order),
                      as.integer(maxorder),
                      as.integer(n),
                      radii=double(n*polyeder$nv),
                      PACKAGE="dti")$radii
    dim(radii) <- c(polyeder$nv,n)
    #
    #   use a sphere of radius level as baseline for the ODF
    #
    radii <- radii+level
    #
    #   to display results in a form that the volumes are comparable,
    #   need volume elements around vertices that have volume proportional to radii,
    #   i.e. a radial extension of radii^(1/3)
    #   odfscale = 1 corresponds to using radii directly for ODF-values
    #   values inbetween are possible
    #
    radii <- radii^(1/odfscale)
    radii <- radii/quantile(apply(radii,2,max),quant)*scale
  }
  if(what %in% c("axis","both")){
    colors <- rainbow(1024,end=2/3)
    ranger <- range(fa)
    ind <- 1024-(fa-ranger[1])/(ranger[2]-ranger[1])*1023
    colorvalues <- rep(colors[ind],rep(2*dim(mix)[1],length(ind)))
    andir <- array(.Fortran("mixandir",
                            as.double(orient),
                            as.double(mix),
                            as.integer(order),
                            as.integer(maxorder),
                            as.integer(n),
                            andir=double(3*n*dim(mix)[1]),
                            PACKAGE="dti")$andir,c(3,dim(mix)[1],n1,n2,n3))
    lcoord <- array(0,c(3,2*dim(mix)[1],n))
    for(i in 1:dim(mix)[1]){
      lcoord[,2*i-1,] <-  andir[,i,]*scale+tmean
      lcoord[,2*i,] <-  -andir[,i,]*scale+tmean
    }          
    dim(lcoord) <- c(3,2*dim(mix)[1]*n)
  }
  dim(tmean) <- c(3,n)
  if(!add) {
    open3d()
    par3d(...)
    rgl.bg(color=bgcolor)
  }
  if(what %in% c("odf","both")) show3dODF(radii,polyeder,centers=tmean,minalpha=minalpha,...)
  if(what %in% c("axis","both"))  rgl.lines(lcoord[1,],lcoord[2,],lcoord[3,],color=colorvalues,size=lwd)
  if(fibers){
    tracks <- tracking(obj,mask=mask,minfa=minfa,maxangle=maxangle)
    dd <- tracks@fibers
    startind <- tracks@startind
    dd <- expandFibers(dd,startind)$fibers
    rgl.lines(dd[,1]+vext[1]/2,dd[,2]+vext[2]/2,dd[,3]+vext[3]/2,
              color=rgb(abs(dd[,4]),abs(dd[,5]),abs(dd[,6])),
              size=lwd)
  }
  if(box) bbox3d()
  if(is.character(title)) {
    title3d(title,color="white",cex=1.5)
  } else {
    if(title) title3d(switch(what,"odf"="estimated ODF"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(what,"odf"="estimated ODF"),"\n")
  if(obj@hmax>1) paste("smoothed with hmax=",obj@hmax,"\n")
  invisible(rgl.cur())
})
##############

setMethod("show3d","dtiIndices", function(obj, 
                                          index=c("fa","ga"), xind=NULL, yind=NULL, zind=NULL, method=1, minfa=0,
                                          bgcolor="black", add=FALSE, lwd=1, box=FALSE, title=FALSE, ...){
  index <- tolower(index)
  ## check index
  index <- match.arg(index)
  #if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  index <- tolower(index) 
  if(!(index%in%c("fa","ga"))) stop("index should be either 'FA' or 'GA'\n")
  cube <- selectCube(xind,yind,zind,obj@ddim,prod(obj@ddim))
  xind <- cube$xind
  yind <- cube$yind
  zind <- cube$zind
  n <- cube$n
  n1 <- cube$n1
  n2 <- cube$n2
  n3 <- cube$n3
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  obj <- obj[xind,yind,zind,drop=FALSE]
  vext <- obj@voxelext
  ind <- switch(index,"fa"=obj@fa[xind,yind,zind], "ga"=obj@ga[xind,yind,zind])
  ind[ind<minfa] <- 0
  ind <- ind*min(vext)
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  andir <- obj@andir[,xind,yind,zind]
  if(method==1) {
    andir <- abs(andir)
    dim(andir) <- c(3,n1*n2*n3)
  } else {
    ind1 <- andir[1,,]<0
    dim(andir) <- c(3,n1*n2*n3)
    andir[,ind1] <- - andir[,ind1]
    andir[2,] <- (1+andir[2,])/2
    andir[3,] <- (1+andir[3,])/2
  }
  colorvalues <- rgb(andir[1,],andir[2,],andir[3,])
  dim(andir) <- c(3,n1,n2,n3)
  andir <- sweep(obj@andir,2:4,ind,"*")
  lcoord <- array(0,c(3,2,n1,n2,n3))
  lcoord[,1,,,] <-  andir/2+tmean[,,,,drop=FALSE]
  lcoord[,2,,,] <-  -andir/2+tmean[,,,,drop=FALSE]
  dim(lcoord) <- c(3,2*n1*n2*n3)
  colorvalues <- c(rbind(colorvalues,colorvalues))
  if(!add) {
    open3d()
    par3d(...)
    rgl.bg(color=bgcolor)
  }
  rgl.lines(lcoord[1,],lcoord[2,],lcoord[3,],color=colorvalues,size=lwd)
  if(box) bbox3d()
  if(is.character(title)) {
    title3d(title,color="white",cex=1.5)
  } else {
    if(title) title3d("Main directions",color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),"Main directions of diffusion estimated from the tensor model\n\n")
  if(box) bbox3d()
  invisible(rgl.cur())
})

##############

setMethod("show3d","dwiQball", function(obj,
                                        xind=NULL, yind=NULL, zind=NULL, level=0, quant=.8, scale=.4, odfscale=1,
                                        bgcolor="black", add=FALSE, subdivide=3, maxobjects=729, minalpha=1, 
                                        box=FALSE, title=FALSE,...){
  #if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(!exists("icosa0")) data("polyeders", envir = environment())

  ## WORKAROUND to make "polyeders" not global variable (for R CMD check) 
  icosa0 <- icosa0
  icosa1 <- icosa1
  icosa2 <- icosa2
  icosa3 <- icosa3
  icosa4 <- icosa4
  ## END WORKAROUND

  if(subdivide<0||subdivide>4) subdivide <- 3
  cube <- selectCube(xind,yind,zind,obj@ddim,maxobjects)
  xind <- cube$xind
  yind <- cube$yind
  zind <- cube$zind
  n <- cube$n
  n1 <- cube$n1
  n2 <- cube$n2
  n3 <- cube$n3
  if (obj@orientation[1]==1) xind <- min(xind)+max(xind)-xind
  if (obj@orientation[2]==3) yind <- min(yind)+max(yind)-yind
  if (obj@orientation[3]==4) zind <- min(zind)+max(zind)-zind
  if(n==0) stop("Empty cube specified")
  cat(" selected cube specified by \n xind=",min(xind),":",max(xind),
      "\n yind=",min(yind),":",max(yind),
      "\n zind=",min(zind),":",max(zind),"\n")
  obj <- obj[xind,yind,zind]
  vext <- obj@voxelext
  tmean <- array(0,c(3,n1,n2,n3))
  tmean[1,,,] <- xind*vext[1]
  tmean[2,,,] <- outer(rep(1,n1),yind)*vext[2]
  tmean[3,,,] <- outer(rep(1,n1),outer(rep(1,n2),zind))*vext[3]
  dim(tmean) <- c(3,n)
  polyeder <- switch(subdivide+1,icosa0,icosa1,icosa2,icosa3,icosa4)
  if(obj@what=="sqrtODF"){
    sphdesign <- design.spheven(obj@order,polyeder$vertices,obj@lambda)$design
    radii <- 0
    for(i in (0:obj@forder)+1){
      sphcoef <- obj@sphcoef[,i,,,]
      dim(sphcoef) <- c(dim(sphcoef)[1],prod(dim(sphcoef)[-1]))
      radii <- radii+(t(sphdesign)%*%sphcoef)^2
    }
    # integration over fourier frequencies (in radial direction)
  } else {
    sphcoef <- obj@sphcoef
    sphdesign <- design.spheven(obj@order,polyeder$vertices,obj@lambda)$design
    dim(sphcoef) <- c(dim(sphcoef)[1],prod(dim(sphcoef)[-1]))
    cat("dim(sphcoef)",dim(sphcoef),"dim(sphdesign)",dim(sphdesign),"\n")
    radii <- t(sphdesign)%*%sphcoef
    radii <- array(pmax(0,radii),dim(radii))
    mradii <- apply(radii,2,mean)
    radii <- sweep(radii,2,mradii,"/")+level
  }
  #  avoid negative ODF's, otherwise scaling by volume produces
  #  strange results
  #
  #   rescale and use a sphere of radius level as baseline for the ODF
  #
  radii[is.na(radii)] <- 0
  #
  #   to display results in a form that the volumes are comparable,
  #   need volume elements around vertices that have volume proportional to radii,
  #   i.e. a radial extension of radii^(1/3)
  #   odfscale = 1 corresponds to using radii directly for ODF-values
  #   values inbetween are possible
  #
  radii <- radii^(1/odfscale)
  radii <- radii/quantile(apply(radii,2,max),quant)*scale
  if(!add) {
    open3d()
    par3d(...)
    rgl.bg(color=bgcolor)
  }
  show3dODF(radii,polyeder,centers=tmean,minalpha=minalpha,...)
  if(box) bbox3d()
  if(is.character(title)) {
    title3d(title,color="white",cex=1.5)
  } else {
    if(title) title3d(switch(obj@what,"ODF"="ODF","wODF"="Weighted ODF","aODF"="alternative ODF","adc"="ADC (Sph. Harmonics)"),color="white",cex=1.5)
  }
  cat("\n rgl-device",rgl.cur(),switch(tolower(obj@what),"odf"="Estimated orientation density function (Qball)","aodf"="Estimated orientation density function (Qball)","adc"="estimated apparent diffusion coefficients (sperical harmonics","wodf"="Estimated orientation density function (Aganji et.al. 2009)"),"\n")
  invisible(rgl.cur())
})

setMethod("show3d","dwiFiber", function(obj, 
                                        add=FALSE, bgcolor="black", box=FALSE, title=FALSE, lwd=1, ...){
  #if(!require(rgl)) stop("Package rgl needs to be installed for 3D visualization")
  if(!add) {
    open3d()
    par3d(...)
    rgl.bg(color=bgcolor)
  }
  dd <- obj@fibers
  startind <- obj@startind
  dd <- expandFibers(dd,startind)$fibers
  rgl.lines(dd[,1],dd[,2],dd[,3],
            color=rgb(abs(dd[,4]),abs(dd[,5]),abs(dd[,6])),
            size=lwd)
  if(box) bbox3d()
  if(is.character(title)) {
    title3d(title,color="white",cex=1.5)
  } else {
    if(title) title3d("Fiber tracks",color="white",cex=1.5)
  }
  
  invisible(rgl.cur())
})

## argument which are not yet decided have "???"
setMethod( "show3d", "dkiTensor", function( obj,
                                            xind = NULL, yind = NULL, zind = NULL,
                                            method = 1, # ???
                                            minfa = .3, # ???
                                            mask = NULL,
                                            level = 0, # ???
                                            quant = .8, # ??? 
                                            scale = .4,# ???
                                            bgcolor = "black",
                                            add = FALSE, 
                                            subdivide = 2, # ???
                                            maxobjects = 729,
                                            what = c( "KT", "DT"),
                                            minalpha = .25, # ???
                                            normalize = NULL, # ???
                                            box = FALSE, # ???
                                            title = FALSE, # ???
                                            ...){
  
  ## first test the OpenGL visualization capabilities! 
  #if( !require( rgl)) stop("Package rgl needs to be installed for 3D visualization")
  
  what <- match.arg( what)
  
  if (!exists("icosa1")) data("polyeders", envir = environment())

  ## WORKAROUND to make "polyeders" not global variable (for R CMD check) 
  icosa0 <- icosa0
  icosa1 <- icosa1
  icosa2 <- icosa2
  icosa3 <- icosa3
  icosa4 <- icosa4
  ## END WORKAROUND
  
  if ( ( subdivide < 0) || ( subdivide > 4)) subdivide <- 2
  polyeder <- switch( subdivide + 1, 
                      icosa0, 
                      icosa1, 
                      icosa2, 
                      icosa3, 
                      icosa4)
  
  cube <- selectCube(xind,yind,zind,obj@ddim,maxobjects)
  xind <- cube$xind
  yind <- cube$yind
  zind <- cube$zind
  n <- cube$n
  n1 <- cube$n1
  n2 <- cube$n2
  n3 <- cube$n3
  
  if ( obj@orientation[ 1] == 1) xind <- min( xind) + max( xind) - xind
  if ( obj@orientation[ 2] == 3) yind <- min( yind) + max( yind) - yind
  if ( obj@orientation[ 3] == 4) zind <- min( zind) + max( zind) - zind
  
  if ( n1*n2*n3 == 0) stop("Empty cube specified")
  cat(" selected cube specified by \n xind=", min( xind), ":", max( xind),
      "\n yind=", min( yind), ":", max( yind),
      "\n zind=", min( zind), ":", max( zind), "\n")
  obj <- obj[ xind, yind, zind]
  
  mask <- if( is.null( mask)) obj@mask else (obj@mask & mask[ xind, yind, zind])
  
  vext <- obj@voxelext
  
  tmean <- array( 0, c( 3, n1, n2, n3))
  tmean[ 1, , , ] <- xind * vext[ 1]
  tmean[ 2, , , ] <- outer( rep( 1, n1), yind) * vext[ 2]
  tmean[ 3, , , ] <- outer( rep( 1, n1), outer( rep( 1, n2), zind)) * vext[ 3]
  dim( tmean) <- c( 3, n)
  
  objind <- dkiIndices( obj)
  andir <- objind@andir
  dim( andir) <- c( 3, n)
  fa <- objind@fa
  if ( method == 1) {
    andir <- abs( andir)
  } else {
    ind <- ( andir[ 1, ] < 0)
    andir[ , ind] <- - andir[ , ind]
    andir[ 2, ] <- ( 1 + andir[ 2, ]) / 2
    andir[ 3, ] <- ( 1 + andir[ 3, ]) / 2
  }
  colorvalues <- rgb( andir[ 1, ], andir[ 2, ], andir[ 3, ])
  
  D <- obj@D
  dim(D) <- c( 6, n)
  
  xxx <- dkiDesign( polyeder$vertices)
  
  Dapp <- xxx[ , c( 1, 4, 5, 2, 6, 3)] %*% D
  
  if ( what == "KT") {
    
    W <- obj@W
    dim(W) <- c( 15, n)
    
    MD <- apply( D[ c( 1, 4, 6),], 2, mean)^2
    
    radii <- sweep( ( xxx[ , 7:21] %*% W) / Dapp, 2 , MD, "*")
    radii[ radii < 0] <- 0
    
    radii <- radii / 2.5 / max( radii) * min( vext)
    
  } else { ## "DT"
    
    radii <- Dapp
    radii[ radii < 0] <- 0
    
    radii <- radii / 2.5 / max( radii) * min( vext)
    
  }
  
  if ( !add) {
    open3d()
    par3d( ...)
    rgl.bg( color = bgcolor)
  }
  show3dTens( radii, polyeder, centers = tmean, colors = colorvalues, alpha = minalpha + ( 1 - minalpha) * fa, ...)
  
  invisible(rgl.cur())
})




################################################################
#                                                              #
# Section for show3d() functions (misc)                        #
#                                                              #
################################################################

show3dTens <- function(radii, polyeder, centers=NULL, colors=NULL, alpha=1, ...){
  if(is.null(centers)){
    centers <- matrix(0,3,1)
    n <- 1
  } else {
    dcenters <- dim(centers)
    if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
    n <- dcenters[2]
  }
  if(is.null(colors)){
    colors <- heat.colors(1)
  } 
  if(length(colors)!=n){
    nc <- length(colors)
    nnc <- n%/%nc+1
    colors <- rep(colors,nnc)[1:n]
  }
  if(is.null(alpha)){
    alpha <- 1
  } 
  if(length(alpha)!=n){
    nc <- length(alpha)
    nnc <- n%/%nc+1
    alpha <- rep(alpha,nnc)[1:n]
  }
  nv <- polyeder$nv
  ni <- polyeder$ni*3
  colors <- t(matrix(colors,n,ni))
  alpha <- t(matrix(alpha,n,ni))
  vertices <- array(polyeder$vertices,c(3,nv,n))
  indices <- matrix(polyeder$indices,c(ni,n))
  if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
  vertices <- sweep(vertices,2:3,radii,"*")
  vertices <- sweep(vertices,c(1,3),centers,"+")
  dim(vertices) <- c(3,nv*n)
  indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
  rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                color=colors,alpha=alpha,...)
}

#############

show3dData <- function( radii, vertices, centers=NULL, minalpha=1, ...){
  #
  #   use gradients directly
  #
  if(is.null(centers)){
    centers <- matrix(0,3,1)
    n <- 1
  } else {
    dcenters <- dim(centers)
    if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
    n <- dcenters[2]
  }
  maxradii <- apply(radii,2,max)
  alpha <- minalpha+(1-minalpha)*sweep(radii,2,maxradii,"/")
  colors <- rgb(abs(vertices[1,]),abs(vertices[2,]),abs(vertices[3,]))
  if(length(alpha)!=n){
    nc <- length(alpha)
    nnc <- n%/%nc+1
    alpha <- rep(alpha,nnc)[1:n]
  }
  nv <- dim(vertices)[2]
  lines <- array(0,c(2,3,n,nv))
  vertices <- array(vertices,c(3,nv,n))
  colors <- matrix(colors,nv,n)
  ind1 <- rep(1:(n*nv),rep(2,n*nv))
  if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(vertices)[2]*dim(centers)[2]")
  vertices0 <- sweep(vertices,2:3,radii,"*")
  vertices <- sweep(.9*vertices0,c(1,3),centers,"+")
  lines[1,,,] <- aperm(vertices,c(1,3,2))
  vertices <- sweep(vertices0,c(1,3),centers,"+")
  lines[2,,,] <- aperm(vertices,c(1,3,2))
  colors <- array(colors,c(nv,n))
  #   rgl.lines(lines[,1,,],lines[,2,,],lines[,3,,],color=t(colors)[ind1],lwd=1)
  rgl.lines(lines[,1,,],lines[,2,,],lines[,3,,],lwd=1)
  rgl.points(vertices[1,,],vertices[2,,],vertices[3,,],color=colors,size=4)
  vertices <- sweep(-.9*vertices0,c(1,3),centers,"+")
  lines[1,,,] <- aperm(vertices,c(1,3,2))
  vertices <- sweep(-vertices0,c(1,3),centers,"+")
  lines[2,,,] <- aperm(vertices,c(1,3,2))
  #   rgl.lines(lines[,1,,],lines[,2,,],lines[,3,,],color=t(colors)[ind1],lwd=1)
  rgl.lines(lines[,1,,],lines[,2,,],lines[,3,,],lwd=1)
  rgl.points(vertices[1,,],vertices[2,,],vertices[3,,],color=colors,size=4)
}

#############

show3dCdata <- function( radii, polyeder, centers=NULL, minalpha=1, scale=.5, ...){
  if(is.null(centers)){
    centers <- matrix(0,3,1)
    n <- 1
  } else {
    dcenters <- dim(centers)
    if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
    n <- dcenters[2]
  }
  vertices <- polyeder$vertices
  colors <- rainbow(1000,start=0,end=.7)[pmax(1,pmin(1000,as.integer(1000*radii)))]
  alpha <- minalpha+(1-minalpha)*radii
  nv <- polyeder$nv
  ni <- polyeder$ni*3
  colors <- matrix(colors,nv,n)
  alpha <- matrix(alpha,nv,n)
  vertices <- array(polyeder$vertices,c(3,nv,n))
  indices <- matrix(polyeder$indices,c(ni,n))
  if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
  vertices <- sweep(scale*vertices,c(1,3),centers,"+")
  dim(vertices) <- c(3,nv*n)
  indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
  rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                color=colors[indices],alpha=alpha[indices],...)
}

#############

show3dODF <- function( radii, polyeder, centers=NULL, minalpha=1, ...){
  if(is.null(centers)||length(centers)==3){
    centers <- matrix(0,3,1)
    n <- 1
  } else {
    dcenters <- dim(centers)
    if(length(dcenters)!=2 || dcenters[1]!=3) stop("centers needs to be NULL or a matrix 
      with dimension (3,n)")
    n <- dcenters[2]
  }
  vertices <- polyeder$vertices
  alpha <- rep(minalpha,length(radii))
  nv <- polyeder$nv
  ni <- polyeder$ni*3
  colors <- rainbow(1024,end=2/3)
  rradii <- apply(radii,2,range)
  sradii <- sweep(sweep(radii,2,rradii[1,],"-"),2,rradii[2,]-rradii[1,],"/")
  ind <- 1024-sradii*1023
  alpha <- matrix(alpha,nv,n)
  vertices <- array(polyeder$vertices,c(3,nv,n))
  indices <- array(polyeder$indices,c(ni,n))
  if(length(radii)!=nv*n) stop("wrong length of radii, needs to be 
             dim(polyeder$vertices)[2]*dim(centers)[2]")
  vertices <- sweep(vertices,2:3,radii,"*")
  vertices <- sweep(vertices,c(1,3),centers,"+")
  dim(vertices) <- c(3,nv*n)
  indices <- sweep(matrix(indices,ni,n),2,((1:n)-1)*nv,"+")
  rgl.triangles(vertices[1,indices],vertices[2,indices],vertices[3,indices],
                color=colors[ind][indices],alpha=alpha[indices],...)
}
