# This file contains the implementation of dti.smooth() for 
# "dtiData" lines 12--
# and 
# "dtiTensor" lines 225--
# objects. The latter can be done with 
# "Euclidian Metric" lines 237--
# "Riemann Metric" lines 414--

dti.smooth <- function(object, ...) cat("No DTI smoothing defined for this class:",class(object),"\n")

setGeneric("dti.smooth", function(object, ...) standardGeneric("dti.smooth"))

setMethod("dti.smooth", "dtiData", function(object,hmax=5,hinit=NULL,lambda=20,tau=10,
                                            rho=1,graph=FALSE,slice=NULL,quant=.8,
                                            minfa=NULL,hsig=2.5,lseq=NULL, method="nonlinear",rician=TRUE,niter=5,result="Tensor") {
switch(method,"linear" = dtilin.smooth(object,hmax,hinit,lambda,rho,graph,slice,quant,minfa,hsig,lseq),
              "nonlinear" =  dtireg.smooth(object,hmax,hinit,lambda,rho,graph,slice,quant,minfa,hsig,lseq,rician,niter,result))
}
)
dtilin.smooth <- function(object,hmax=5,hinit=NULL,lambda=52,
                                            rho=1,graph=FALSE,slice=NULL,quant=.8,
                                            minfa=NULL,hsig=2.5,lseq=NULL){
#
#     lambda and lseq adjusted for alpha=0.2
#
  wlse <- TRUE
  eps <- 1e-6
#  if (graph) {
#    adimpro <- require(adimpro)
#    if (!adimpro) cat("No graphical output! Install package adimpro from CRAN!\n")
#    graph <- graph & adimpro
#  }
  args <- sys.call(-3)
  args <- c(object@call,args)
  s0ind <- object@s0ind
  ns0 <- length(s0ind)
  ngrad <- object@ngrad
  ddim0 <- object@ddim0
  ddim <- object@ddim
  nvox <- prod(ddim)
  z <- sioutlier(object@si,s0ind)
  si <- aperm(array(z$si,c(ngrad,ddim)),c(2:4,1))
  object@si <- si
  index <- z$index
  s0 <- si[,,,s0ind]
  if(length(s0ind)>1) s0 <- apply(s0,1:3,mean) 
  si <- si[,,,-s0ind]
  xind <- object@xind
  yind <- object@yind
  zind <- object@zind
  source <- object@source
  btb <- object@btb[,-s0ind]
  voxelext <- object@voxelext
  if(is.null(voxelext)) vext <- c(1,1,1) else vext <- voxelext/min(voxelext)
  Bcov <- btb%*%t(btb)
  btbsvd <- svd(btb)
  solvebtb <- btbsvd$u %*% diag(1/btbsvd$d) %*% t(btbsvd$v)
  dtobject <- dtiTensor(object,method="linear")
  scale <- dtobject@scale
  mask <- dtobject@mask
  y <- dtobject@D
  sigma2 <- dtobject@sigma
  scorr <- dtobject@scorr
  h0 <- dtobject@bw
  th0 <- dtobject@th0
  cat("Corresponding bandwiths for specified correlation:",h0,"\n")

  rm(object,dtobject)
  gc()
  dimy <- dim(y)
  if(length(dimy)!=4||dimy[1]!=6) stop("y does not contain 3D diffusion tensor image")
  n1<-dimy[2]
  n2<-dimy[3]
  n3<-dimy[4]
  n<-n1*n2*n3
  sigma2[sigma2<=mean(sigma2)*1e-5]<- mean(sigma2)*1e-5
  z <- .Fortran("projdt2",
                as.double(y),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                D=double(6*n),
                anindex=double(n),
                andirection=double(3*n),
                det=double(n),
                as.double(eps),
                PACKAGE="dti")[c("D","anindex","andirection","det")]
  y <- array(z$D,dimy)
  z$bi <- 1/sigma2
  dim(z$D) <- dimy
  dim(z$anindex) <-dim(z$det) <- dimy[-1]
  dim(z$andirection) <- c(3,dimy[-1]) 
  sigma2hat <- .Fortran("smsigma",
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.double(hsig),
                       as.double(vext),
                       sigma2hat=double(n1*n2*n3),
                       PACKAGE="dti")$sigma2hat
   z$sigma2hat <- sigma2hat
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
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxy: min",signif(min(z$D[2,,,][mask]),3),"max",signif(max(z$D[2,,,][mask]),3)))
     img<-z$D[3,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
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
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dyz: min",signif(min(z$D[5,,,][mask]),3),"max",signif(max(z$D[5,,,][mask]),3)))
     andir2.image(z,slice,quant=quant,minfa=minfa,xaxt="n",yaxt="n")
     title(paste("Directions (h=1), slice",slice))
     ni<-array(1,dimy[-1])*as.integer(mask)
     show.image(make.image((65535*ni/max(ni))[,,slice]),xaxt="n",yaxt="n")
         title(paste("sum of weights  mean=",signif(1,3)))
     img<-z$D[6,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dzz: mean",signif(mean(z$D[6,,,][mask]),3),"max",signif(max(z$D[6,,,][mask]),3)))
  }
  ngrad <- ngrad-length(s0ind)
  hincr <- 1.25^(1/3)
  if(is.null(hinit)){
  hakt0 <- 1
  hakt <- hincr
  hinit <- 1
  } else {
  hakt0 <- max(1,hinit/hincr)
  hakt <- hinit
  }
  steps <- as.integer(log(hmax/hinit)/log(hincr)+1)

  # define lseq
  if (is.null(lseq)) {
# this is optimized for lkern="Gaussian" such that alpha approx 0.04 -- 0.1 and probability of separated points is approx. 1e-4
    lseq <- c(1.5,.9,.8,.85,.8,0.85,.85,0.95,1.25,1.15)# alpha=0.2    abs. deviation in S0
  }
  if (length(lseq)<steps) lseq <- c(lseq,rep(1,steps-length(lseq)))
  lseq <- lseq[1:steps]
  k <- 1
  lambda0 <- lambda*lseq[k]
  while(hakt <= hmax) {
    if (any(h0 >= 0.25)) {
       corrfactor <- Spatialvar.gauss(hakt0/0.42445/4,h0,3) /
       Spatialvar.gauss(h0,1e-5,3) /
       Spatialvar.gauss(hakt0/0.42445/4,1e-5,3)
       lambda0 <- lambda0 * corrfactor
       cat("Correction factor for spatial correlation",signif(corrfactor,3),"\n")
    }
     z <- .Fortran("awssidti",
                as.double(s0),
                as.double(si),
                as.logical(mask),
                as.double(z$D),
                bi=as.double(z$bi),
                anindex=as.double(z$anindex),
                andirection=as.double(z$andirection),
                det=as.double(z$det),
                as.double(Bcov),
                as.double(solvebtb),
                as.double(sigma2),
                as.double(z$sigma2hat),
                as.integer(n1),
                as.integer(n2),
                as.integer(n3),
                as.integer(ngrad),
                as.double(hakt),
                as.double(vext),
                as.double(rho),
                as.double(lambda0),
                D=double(6*n),
                sigma2hat=double(n),
                double(ngrad),
                as.double(eps),
                PACKAGE="dti")[c("D","bi","anindex","andirection","det","sigma2hat")]
     if(hakt<hsig){
        eta <- (hsig^3 - hakt^3)/hsig^3
        z$sigma2hat <- eta*sigma2hat+(1-eta)*z$sigma2hat
     }
     dim(z$bi) <- dim(z$anindex) <- dim(z$det) <- dim(z$sigma2hat) <- dimy[-1]
     dim(z$D) <- dimy
     dim(z$andirection) <- c(3,dimy[-1]) 
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
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dxy: min",signif(min(z$D[2,,,][mask]),3),"max",signif(max(z$D[2,,,][mask]),3)))
     img<-z$D[3,,,slice]
     rg<-quantile(img,c(.01,.99))
     img[img>rg[2]]<-rg[2]
     img[img<rg[1]]<-rg[1]
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
     show.image(make.image(img),xaxt="n",yaxt="n")
     title(paste("Dyz: min",signif(min(z$D[5,,,][mask]),3),"max",signif(max(z$D[5,,,][mask]),3)))
     andir2.image(z,slice,quant=quant,minfa=minfa,xaxt="n",yaxt="n")
     title(paste("Directions (h=",signif(hakt,3),"), slice",slice))
     ni<-(z$bi[,,slice]*if(wlse) z$sigma2hat[,,slice] else 1)*mask[,,slice]
     show.image(make.image(65535*ni/max(ni)),xaxt="n",yaxt="n")
     title(paste("sum of weights  mean=",signif(mean((z$bi*if(wlse) z$sigma2hat else 1)[mask]),3)))
     img<-z$D[6,,,slice]
     show.image(make.image(65535*img/max(img)),xaxt="n",yaxt="n")
     title(paste("Dzz: mean",signif(mean(z$D[6,,,][mask]),3),"max",signif(max(z$D[6,,,][mask]),3)))
     }
     cat("h=",signif(hakt,3),"Quantiles (.5, .75, .9, .95, 1) of anisotropy index",signif(quantile(z$anindex,c(.5, .75, .9, .95, 1)),3),"\n")
     hakt0<-hakt
     hakt <- hakt*hincr
    c1 <- (prod(h0+1))^(1/3)
    c1 <- 2.7214286 - 3.9476190*c1 + 1.6928571*c1*c1 - 0.1666667*c1*c1*c1
    x <- (prod(1.25^(k-1)/c(1,1,1)))^(1/3)
    scorrfactor <- (c1+x)/(c1*prod(h0+1)+x)
    k <- k+1
    lambda0 <- lambda*lseq[k]*scorrfactor     
  }
#  replace non-tensors (with negative eigenvalues) by a small isotropic tensor 
      ind <- array(dti3Dev(z$D,mask),c(3,n1,n2,n3))[1,,,]<1e-6
       if(sum(ind&mask)>0){
           dim(z$D) <- c(6,n1*n2*n3)
           z$D[c(1,4,6),ind&mask] <- 1e-6
           z$D[c(2,3,5),ind&mask] <- 0
           dim(z$D) <- c(6,n1,n2,n3)
       }
  
  invisible(new("dtiTensor",
                list( hmax = hmax, ni=z$bi),
                call = args,
                D     = z$D,
                th0   = th0,
                gradient = object@gradient,
                bvalue = object@bvalue,
                btb   = btb,
                sigma = z$sigma2hat, 
                scorr = scorr,
                bw    = h0, 
                mask  = mask,
                hmax  = 1,
                ngrad = ngrad+length(s0ind), # = dim(btb)[2]
                s0ind = s0ind,
                ddim  = as.integer(ddim),
                ddim0 = as.integer(ddim0),
                xind  = xind,
                yind  = yind,
                zind  = zind,
                voxelext = voxelext,
                orientation = object@orientation,
                rotation = object@rotation,
                source= source,
                outlier = index,
                scale = scale,
                method= "linear")
            )

}

