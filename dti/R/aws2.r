#
#   Nonlinear regression; regularized version according to Koay et. al. (2006)
#   this is also based on an a statistical penalty defined using log-likelihood difference
#
dtireg.smooth <- function(object,hmax=5,hinit=1,lambda=30,rho=1,graph=FALSE,slice=NULL,quant=.8,
                         minfa=NULL,hsig=2.5,lseq=NULL,rician=TRUE,niter=5,result="Tensor"){
#
#     lambda and lseq adjusted for alpha=0.2
#
  if(!is.null(object$ni)){
     warning("DWI object has been smoothed already, smoothing omitted")
     return(if(result=="Tensor") dtiTensor(object,method="nonlinear") else object)
  }
  eps <- 1e-6
  maxnw <- 10000
#  if (graph) {
#    adimpro <- require(adimpro)
#    if (!adimpro) cat("No graphical output! Install package adimpro from CRAN!\n")
#    graph <- graph & adimpro
#  }
  args <- sys.call(-3)
  args <- c(object@call,args)
  sdcoef <- object@sdcoef
  dtobject <- dtiTensor(object,method="nonlinear")
  scale <- dtobject@scale
  mask <- dtobject@mask
  th0 <- dtobject@th0
  D <- dtobject@D
  scorr <- dtobject@scorr
  h0 <- dtobject@bw
  rm(dtobject)
  gc()
  cat("Corresponding bandwiths for specified correlation:",h0,"\n")
  s0ind <- object@s0ind
  ngrad <- object@ngrad
  ddim0 <- object@ddim0
  ddim <- object@ddim
  z <- sioutlier(object@si,s0ind)
  si <- array(z$si,c(ngrad,ddim))
  index <- z$index 
  rm(z)
  gc()
  xind <- object@xind
  yind <- object@yind
  zind <- object@zind
  source <- object@source
  btb <- object@btb
  voxelext <- object@voxelext
  if(is.null(voxelext)) vext <- c(1,1,1) else vext <- voxelext/min(voxelext)
  gc()
  dimy <- dim(D)
  if(length(dimy)!=4||dimy[1]!=6) stop("D does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  z <- .Fortran("projdt2",
                as.double(D),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                D=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                as.double(eps),
                PACKAGE="dti")[c("D","anindex","andirection","det")]
  dim(z$D) <- dimy
  z$th0 <- th0
  dim(z$anindex) <-dim(z$det) <- dimy[-1]
  dim(z$andirection) <- c(3,dimy[-1]) 
  z$sihat <- si
  z$bi <- rep(1,n)
#
#  initial state for h=1
#
  if(graph){
     oldpar <- par(mfrow=c(3,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
     on.exit(par(oldpar))
     if(is.null(slice)) slice<-n3%/%2
     class(z) <- "dti"
     img<-z$D[1,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dxx: mean",signif(mean(z$D[1,,,][mask]),3),"max",signif(max(z$D[1,,,][mask]),3)))
     img<-z$D[2,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     img <- (img-rg[1])/(rg[2]-rg[1])
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxy: min",signif(min(z$D[2,,,][mask]),3),"max",signif(max(z$D[2,,,][mask]),3)))
     img<-z$D[3,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     img <- (img-rg[1])/(rg[2]-rg[1])
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxz: min",signif(min(z$D[3,,,][mask]),3),"max",signif(max(z$D[3,,,][mask]),3)))
     show.image(make.image(z$anindex[,,slice]),xaxt="n",yaxt="n")
     title(paste("Anisotropy index  range:",signif(min(z$anindex[mask]),3),"-",
                  signif(max(z$anindex[mask]),3)))
     img<-z$D[4,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dyy: mean",signif(mean(z$D[4,,,][mask]),3),"max",signif(max(z$D[4,,,][mask]),3)))
     img<-z$D[5,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     img <- (img-rg[1])/(rg[2]-rg[1])
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dyz: min",signif(min(z$D[5,,,][mask]),3),"max",signif(max(z$D[5,,,][mask]),3)))
     andir2.image(z,slice,quant=quant,minfa=minfa,xaxt="n",yaxt="n")
     title(paste("Directions (h=1), slice",slice))
     ni <- array(1,dimy[-1])*as.integer(mask)
     show.image(make.image((65535*ni/max(ni))[,,slice]),xaxt="n",yaxt="n")
     title(paste("sum of weights  mean=",signif(1,3)))
     img<-z$D[6,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dzz: mean",signif(mean(z$D[6,,,][mask]),3),"max",signif(max(z$D[6,,,][mask]),3)))
  }
  hincr <- 1.25^(1/3)
  maxvol <- getvofh(hmax,c(1,0,0,1,0,1),vext)
  kstar <- as.integer(log(maxvol)/log(1.25)) 
  if(is.null(hinit)){
  hakt0 <- 1
  hakt <- hincr
  hinit <- 1
  } else {
  hakt0 <- max(1,hinit/hincr)
  hakt <- hinit
  }
  steps <- kstar+1

  # define lseq
  if (is.null(lseq)) {
# this is optimized for lkern="Gaussian" such that alpha approx 0.04 -- 0.1 and probability of separated points is approx. 1e-4
    lseq <- c(1.5,.9,.8,.85,.8,0.85,.85,0.95,1.25,1.15)# alpha=0.2    abs. deviation in S0
  }
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  k <- 1
  lambda0 <- lambda
  while(k <= kstar) {
      hakt0 <- gethani(1,10,1.25^(k-1),c(1,0,0,1,0,1),vext,1e-4)
      hakt <- gethani(1,10,1.25^k,c(1,0,0,1,0,1),vext,1e-4)
    if (any(h0 >= 0.25)) {
       corrfactor <- Spatialvar.gauss(hakt0/0.42445/4,h0,3) /
       Spatialvar.gauss(h0,1e-5,3) /
       Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
       lambda0 <- lambda0 * corrfactor
       cat("Correction factor for spatial correlation",signif(corrfactor,3),"\n")
    }
    z <- .Fortran("awsrgdti",
                    as.double(si),
                    sihat=as.double(z$sihat), # needed for statistical penalty
                    double(ngrad*n),# array for predicted Si's from the tensor model 
                    as.integer(ngrad),
                    as.integer(n1),
                    as.integer(n2),
                    as.integer(n3),
                    as.logical(mask),
                    as.logical(!((1:ngrad)%in%s0ind)),
                    as.double(btb),
                    as.double(sdcoef),
                    as.double(z$th0),
                    th0=double(n),
                    as.double(z$D),
                    D=double(6*n),
                    bi=as.double(z$bi),
                    anindex=as.double(z$anindex),
                    andirection=as.double(z$andirection),
                    det=as.double(z$det),
                    sigma2r=double(n),
                    as.double(1.25^k),
                    as.integer(niter),
                    as.double(vext),
                    as.double(rho),
                    as.double(lambda0),
                    double(ngrad),#swsi
                    double(ngrad),#swsi2
                    double(ngrad),#F
                    double(ngrad),#var
                    as.double(eps),
                    as.logical(rician), 
                    as.integer(maxnw),# maximum number of positive weights
                    integer(ngrad),# auxiliary for number of iterations
                    double(maxnw*ngrad),# auxiliary for aktive data
                    integer(maxnw*3),# auxiliary for index of aktive data
                    double(maxnw),# auxiliary for weights
                    double(ngrad),# auxiliary for variances
                    PACKAGE="dti")[c("th0","D","bi","anindex","andirection","det","sigma2r","sihat")] 
     dim(z$th0) <- dim(z$bi) <- dim(z$anindex) <- dim(z$det) <- dim(z$sigma2r) <- dimy[-1]
     dim(z$D) <- dimy
     dim(z$andirection) <- c(3,dimy[-1]) 
     n <- n1*n2*n3
     if(any(is.na(z$th0))){
        indna <- is.na(z$th0)
#        cat("found ",sum(indna),"NA's\n")
#        cat("indna",(1:n)[indna],"\n")
        z$th0[indna] <- th0[indna]
#        cat("th0",th0[indna],"\n")
        dim(z$D) <- dim(D) <- c(6,n1*n2*n3)
        z$D[,indna] <- D[,indna]
        dim(z$D) <- dim(D) <- c(6,n1,n2,n3)
        mask[indna] <- FALSE
     }
##
##  set mask to FALSE on voxel where we observe extreme 
##  values of z$det
##  this is an indication for numerical problems caused
##  e.g. by registration artifacts
##
     det95 <- quantile(z$det[mask],.95)
     if(any(z$det[mask]>1e3*det95)){
        indna <- (1:n)[z$det/1e3>det95&mask]
#        cat("found ",length(indna),"voxel with extreme derterminat 's, keep initial estimates and do not use these voxel\n")
#        cat("indna",indna,"\n")
        z$th0[indna] <- th0[indna]
        dim(z$D) <- dim(D) <- c(6,n1*n2*n3)
        z$D[,indna] <- D[,indna]
        dim(z$D) <- dim(D) <- c(6,n1,n2,n3)
        dim(z$sihat)  <- dim(si) <- c(ngrad,n1*n2*n3)
        z$sihat[,indna] <- si[,indna]
        dim(z$sihat) <- dim(si) <- c(ngrad,n1,n2,n3)
        z$det[indna] <- 0
        mask[indna] <- FALSE
     }
     if(graph){
     class(z) <- "dti"
     img<-z$D[1,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dxx: mean",signif(mean(z$D[1,,,][mask]),3),"max",signif(max(z$D[1,,,][mask]),3)))
     img<-z$D[2,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     img <- (img-rg[1])/(rg[2]-rg[1])
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxy: min",signif(min(z$D[2,,,][mask]),3),"max",signif(max(z$D[2,,,][mask]),3)))
     img<-z$D[3,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     img <- (img-rg[1])/(rg[2]-rg[1])
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxz: min",signif(min(z$D[3,,,][mask]),3),"max",signif(max(z$D[3,,,][mask]),3)))
     show.image(make.image(z$anindex[,,slice]),xaxt="n",yaxt="n")
     title(paste("Anisotropy index  range:",signif(min(z$anindex[mask]),3),"-",
                  signif(max(z$anindex[mask]),3)))
     img<-z$D[4,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dyy: mean",signif(mean(z$D[4,,,][mask]),3),"max",signif(max(z$D[4,,,][mask]),3)))
     img<-z$D[5,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
     img <- (img-rg[1])/(rg[2]-rg[1])
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dyz: min",signif(min(z$D[5,,,][mask]),3),"max",signif(max(z$D[5,,,][mask]),3)))
     andir2.image(z,slice,quant=quant,minfa=minfa,xaxt="n",yaxt="n")
     title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
     ni<-z$bi[,,slice]*mask[,,slice]
     show.image(make.image(65535*ni/max(ni)),xaxt="n",yaxt="n")
     title(paste("sum of weights  mean=",signif(mean(z$bi[mask]),3)))
     img<-z$D[6,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dzz: mean",signif(mean(z$D[6,,,][mask]),3),"max",signif(max(z$D[6,,,][mask]),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex[mask],c(.5, .75, .9, .95, 1)),3),"\n")
    k <- k+1
     lambda0 <- lambda 
     gc()
  }
  dimsi <- dim(si)
  rm(si)
  gc()
  if(result=="Tensor"){
  cat("prepare final dtiTensor object",format(Sys.time()),"\n")
  } else   cat("prepare final smoothed dtiData object",format(Sys.time()),"\n")
  if(result=="Tensor") invisible(new("dtiTensor",
                list(s2rician=if(rician) z$sigma2r else NULL, ni=z$bi),
                call = args,
                D = z$D,
                th0 = z$th0,
                sigma = array(0,c(1,1,1)),
                scorr = scorr,
                bw = h0,
                mask = mask,
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb   = btb,
                hmax  = hmax,
                ngrad = ngrad, # = dim(btb)[2]
                s0ind = object@s0ind,
                ddim  = as.integer(ddim),
                ddim0 = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                voxelext = object@voxelext,
                orientation = object@orientation,
                rotation = object@rotation,
                source= object@source,
                outlier = index,
                scale = scale,
                method= "nonlinear")
            ) else invisible(new("dtiData",
                list(s2rician=if(rician) z$sigma2r else NULL, ni=z$bi),
                call = args,
                si = aperm(array(z$sihat,dimsi),c(2:4,1)),
                sdcoef = sdcoef,
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb    = btb,
                ngrad  = ngrad, # = dim(btb)[2]
                s0ind  = object@s0ind, # indices of s0 images
                replind = object@replind, # replications in gradient design
                ddim   = as.integer(ddim),
                ddim0  = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                voxelext = object@voxelext,
                level  = object@level,
                orientation = object@orientation,
                rotation = object@rotation,
                source= object@source)
            )
}

gethani <- function(x,y,value,a,vext,eps=1e-2){
.Fortran("gethani",
         as.double(x),
         as.double(y),
         as.double(value),
         as.double(a),
         as.double(vext),
         as.double(eps),
         bw=double(1),
         PACKAGE="dti")$bw
}
getvofh <- function(bw,a,vext){
.Fortran("getvofh",
         as.double(a),
         as.double(bw),
         as.double(vext),
         vol=double(1),
         PACKAGE="dti")$vol
}

