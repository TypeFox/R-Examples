#########################################################################################################################
#
#    R - function  aws3D  for vector-valued  Adaptive Weights Smoothing (AWS) in 3D
#    
#    reduced version for qtau==1, heteroskedastic variances only, exact value for Variance reduction
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2005 Weierstrass-Institut fuer
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  Joerg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#

vaws3D <- function(y,qlambda=NULL,lkern="Gaussian",skern="Plateau",weighted=TRUE,
                   sigma2=NULL,mask=NULL,hinit=NULL,hincr=NULL,hmax=NULL,
                   ladjust=1,u=NULL,graph=FALSE,demo=FALSE,wghts=NULL,
                   spmin=.3,h0=c(0,0,0),vwghts=1,testprop=FALSE,
                   res=NULL, resscale=NULL, ddim=NULL){
#
#  qlambda, corrfactor adjusted for case lkern="Gaussian",skern="Plateau" only
#
#  uses sum of weights and correction factor C(h,g) in statistical penalty
#
  #  Auxilary functions
  IQRdiff <- function(y) IQR(diff(y))/1.908

  #
  # first check arguments and initialize
  #
  # test dimension of data (vector of 3D) and define dimension related stuff
  if(is.null(res)) {
      return(warning("Please specify keep=''all'' when calling fmri.lm"))
  }
  d <- 3
  dy <- dim(y)
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1*n2*n3
  if (length(dy)==d) {
    dim(y) <- dy <- c(dy,1)
  } else if (length(dy)!=d+1) {
    stop("y has to be 3 or 4 dimensional")
  }
  dv <- dim(y)[d+1]
  # define vwghts
  if (length(vwghts)>dv) vwghts <- vwghts[1:dv]
  dv0 <- sum(vwghts>0)

  # MAE
  mae <- NULL

  # set the code for the kernel (used in lkern) and set lambda
  lkern <- switch(lkern,
                  Triangle=2,
                  Plateau=1,
                  Gaussian=3,
                  1)
  skern <- switch(skern,
                  Triangle=2,
                  Plateau=1,
                  Exponential=3,
                  2)

  # define qlambda, lambda
  if (is.null(qlambda)) qlambda <- switch(skern,.9945,.9988,.9995)
  if (qlambda<.9) warning("Inappropriate value of qlambda")
  if (qlambda<1) {
    lambda <- ladjust*qchisq(qlambda,dv0) 
  } else {
    lambda <- 1e50
  }
  if(skern%in%c(1,2)) {
    # to have similar preformance compared to skern="Exp"
    lambda <- 4/3*lambda
    if(skern==1) if(is.null(spmin)) spmin <- .3 
    spmax <- 1
  } else {
      spmin <- 0
      spmax <- 4
  }
  # set hinit and hincr if not provided
  if (is.null(hinit)||hinit<1) hinit <- 1
  
  # define hmax
  if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points

  # re-define bandwidth for Gaussian lkern!!!!
  if (lkern==3) {
    # assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- fwhm2bw(hmax)*4
    hinit <- min(hinit,hmax)
  }
  if (qlambda == 1) hinit <- hmax

  # define hincr
  if (is.null(hincr)||hincr<=1) hincr <- 1.25
  hincr <- hincr^(1/d)

  # determine corresponding bandwidth for specified correlation
  if(is.null(h0)) h0 <- rep(0,3)

  # estimate variance in the gaussian case if necessary  
  if (is.null(sigma2)) {
    sigma2 <- IQRdiff(as.vector(y))^2
    if (any(h0>0)) sigma2<-sigma2*Varcor.gauss(h0)*Spatialvar.gauss(h0,1e-5,d)
    if (demo) cat("Estimated variance: ", signif(sigma2,4),"\n")
  }
  if (length(sigma2)==1) sigma2<-array(sigma2,dy[1:3]) 
  if(is.null(mask)) mask <- array(TRUE,dy[1:3])
  mask[sigma2>=1e16] <- FALSE
#  in these points sigma2 probably contains NA's
  # deal with homoskedastic Gaussian case by extending sigma2
  if (length(sigma2)!=n) stop("sigma2 does not have length 1 or same length as y")
  dim(sigma2) <- dy[1:3]
  lambda <- lambda*2 
  sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 

  # demo mode should deliver graphics!
  if (demo && !graph) graph <- TRUE
  
  # Initialize  list for bi and theta
  if (is.null(wghts)) wghts <- c(1,1,1)
  hinit <- hinit/wghts[1]
  hmax <- hmax/wghts[1]
  wghts <- (wghts[2:3]/wghts[1])
  tobj <- list(bi= rep(1,n))
  theta <- y
  maxvol <- getvofh(hmax,lkern,wghts)
  kstar <- as.integer(log(maxvol)/log(1.25))
  steps <- kstar+1
  
  if(qlambda < 1) k <- 1 else k <- kstar
  hakt <- hinit
  hakt0 <- hinit
  lambda0 <- lambda 
  if (hinit>1) lambda0 <- 1e50 # that removes the stochstic term for the first step
  scorr <- numeric(3)
  if(h0[1]>0) scorr[1] <-  get.corr.gauss(h0[1],2)
  if(h0[2]>0) scorr[2] <-  get.corr.gauss(h0[2],2)
  if(h0[3]>0) scorr[3] <-  get.corr.gauss(h0[3],2)
  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  if(testprop) {
    if(is.null(u)) u <- 0
  } 
  propagation <- NULL
  # run single steps to display intermediate results
  while (k<=kstar) {
      hakt0 <- gethani(1,3,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,3,lkern,1.25^k,wghts,1e-4)
      hakt.oscale <- if(lkern==3) bw2fwhm(hakt/4) else hakt
      cat("step",k,"bandwidth",signif(hakt.oscale,3)," ")
    dlw <- (2*trunc(hakt/c(1,wghts))+1)[1:d]
#  need bandwidth in voxel for Spaialvar.gauss, h0 is in voxel
    if (any(h0>0)) lambda0 <- lambda0 * Spatialvar.gauss(bw2fwhm(hakt0)/4/c(1,wghts),h0,d)/
      Spatialvar.gauss(h0,1e-5,d)/Spatialvar.gauss(bw2fwhm(hakt0)/4/c(1,wghts),1e-5,d)
        # Correction C(h0,hakt) for spatial correlation depends on h^{(k-1)}  all bandwidth-arguments in FWHM 
    hakt0 <- hakt
    theta0 <- theta
    bi0 <- tobj$bi
    #
    #   need these values to compute variances after the last iteration
    #
    tobj <- .Fortran("chaws2",as.double(y),
                     as.double(sigma2),
                     as.logical(!mask),
                     as.logical(weighted),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.integer(dv),
                     as.integer(dv0),
                     hakt=as.double(hakt),
                     as.double(lambda0),
                     as.double(theta0),
                     bi=as.double(bi0),
                     thnew=double(n1*n2*n3*dv),
                     as.integer(lkern),
                     as.integer(skern),
                     as.double(spmin),
                     as.double(spmax),
                     double(prod(dlw)),
                     as.double(wghts),
                     as.double(vwghts),
                     double(dv),#swjy
                     double(dv0),#thi
                     PACKAGE="fmri",DUP=TRUE)[c("bi","thnew","hakt")]
    gc()
    theta <- array(tobj$thnew,dy) 
    dim(tobj$bi) <- dy[-4]
    if(testprop) {
      pobj <- .Fortran("chaws2",as.double(y),
                       as.double(sigma2),
                       as.logical(!mask),
                       as.logical(weighted),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
		       as.integer(dv),
		       as.integer(dv0),
                       hakt=as.double(hakt),
                       as.double(1e50),
                       as.double(theta0),
                       bi=as.double(bi0),
                       thnew=double(n1*n2*n3*dv),
                       as.integer(lkern),
		       as.integer(skern),
	               as.double(spmin),
		       as.double(spmax),
		       double(prod(dlw)),
		       as.double(wghts),
		       as.double(vwghts),
		       double(dv),#swjy
		       double(dv0),#thi
		       PACKAGE="fmri",DUP=TRUE)[c("bi","thnew")]
      ptheta <- array(pobj$thnew,dy) 
      rm(pobj) 
      gc()
      propagation <- c(propagation,sum(abs(theta-ptheta))/sum(abs(ptheta-u)))
      cat("Propagation with alpha=",max(propagation),"\n")
      cat("alpha values:","\n")
      print(rbind(hakt.oscale,signif(propagation[-1],3)))
    }
    if (graph) {
      par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
      image(y[,,n3%/%2+1,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
      title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
      image(theta[,,n3%/%2+1,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
      title(paste("Reconstruction  h=",signif(hakt.oscale,3)," min=",signif(min(theta),3)," max=",signif(max(theta),3)))
      image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
      title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
    }
    if (!is.null(u)) {
      cat("bandwidth: ",signif(hakt.oscale,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
          signif(mean((theta-u)^2),3),"   MAE: ",signif(mean(abs(theta-u)),3)," mean(bi)=",signif(mean(tobj$bi),3),"\n")
      mae <- c(mae,signif(mean(abs(theta-u)),3))
    } else if (max(total) >0) {
      cat(signif(total[k],2)*100,"%                 \r",sep="")
     }
    if (demo) readline("Press return")
    k <- k+1
#  adjust lambda for the high intrinsic correlation between  neighboring estimates 
    c1 <- (prod(h0+1))^(1/3)
    c1 <- 2.7214286 - 3.9476190*c1 + 1.6928571*c1*c1 - 0.1666667*c1*c1*c1
    x <- (prod(1.25^(k-1)/c(1,wghts)))^(1/3)
    scorrfactor <- (c1+x)/(c1*prod(h0+1)+x)
    lambda0 <- lambda*scorrfactor 
    gc()
  }

  #   Now compute variance of theta and variance reduction factor (with respect to the spatially uncorrelated situation   
  residuals <- readBin(res,"integer",prod(ddim),2)
  cat("\nfmri.smooth: first variance estimate","\n")
  vartheta0 <- .Fortran("ivar",as.double(residuals),
                           as.double(resscale),
                           as.logical(mask),
                           as.integer(ddim[1]),
                           as.integer(ddim[2]),
                           as.integer(ddim[3]),
                           as.integer(ddim[4]),
                           var = double(n1*n2*n3),
                           PACKAGE="fmri",DUP=TRUE)$var
  cat("fmri.smooth: smooth the residuals","\n")
  residuals <- .Fortran("ihaws2",as.double(residuals),
                     as.double(sigma2),
                     as.logical(!mask),
                     as.logical(weighted),
                     as.integer(ddim[1]),
                     as.integer(ddim[2]),
                     as.integer(ddim[3]),
                     as.integer(ddim[4]),
                     as.integer(dv0),
                     hakt=as.double(tobj$hakt),
                     as.double(lambda0),
                     as.double(theta0),
                     as.double(bi0),
                     resnew=double(prod(ddim)),
                     as.integer(lkern),
                     as.integer(skern),
                     as.double(spmin),
                     as.double(spmax),
                     double(prod(dlw)),
                     as.double(wghts),
                     as.double(vwghts),
                     double(ddim[4]),#swjy
                     double(dv0),#thi
                     PACKAGE="fmri",DUP=TRUE)$resnew
  dim(residuals) <- c(ddim[4],n1,n2,n3)
  gc()
#   get variances and correlations
  cat("fmri.smooth: estimate correlations","\n")
  lags <- c(5,5,3)
  scorr <- .Fortran("imcorr",as.double(residuals),
                     as.logical(mask),
                     as.integer(n3),
                     as.integer(n2),
                     as.integer(n1),
                     as.integer(ddim[4]),
                     scorr=double(prod(lags)),
                     as.integer(lags[1]),
                     as.integer(lags[2]),
                     as.integer(lags[3]),
                     PACKAGE="fmri",DUP=TRUE)$scorr
  dim(scorr) <- lags                     
  cat("fmri.smooth: final variance estimate","\n")
  vartheta <- .Fortran("ivar",as.double(residuals),
                           as.double(resscale),
                           as.logical(mask),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(n3),
                           as.integer(ddim[4]),
                           var = double(n1*n2*n3),
                           PACKAGE="fmri",DUP=TRUE)$var
  vred <- array(vartheta/vartheta0,c(n1,n2,n3))
  vartheta <- vred/sigma2  #  sigma2 contains inverse variances
  z <- list(theta=theta,ni=tobj$bi,var=vartheta,vred=vred,vred0=median(vred[mask]),y=y,
            hmax=tobj$hakt*switch(lkern,1,1,bw2fwhm(1/4)),mae=mae,
            alpha=propagation,scorr=scorr,res=residuals,mask=mask)
  #   vred accounts for variance reduction with respect to uncorrelated (\check{sigma}^2) data
  class(z) <- "aws.gaussian"
  invisible(z)
}

vaws3Dfull <- function(y,qlambda=NULL,lkern="Gaussian",skern="Plateau",weighted=TRUE,
                   sigma2=NULL,mask=NULL,hinit=NULL,hincr=NULL,hmax=NULL,
                   ladjust=1,u=NULL,graph=FALSE,demo=FALSE,wghts=NULL,
                   spmin=.3,vwghts=1,testprop=FALSE,
                   res=NULL, resscale=NULL, ddim=NULL) {
#
#  qlambda, corrfactor adjusted for case lkern="Gaussian",skern="Plateau" only
#
#  implements estimation of true variances from residuals (slower) than vaws3D
#
  #  Auxilary functions
  IQRdiff <- function(y) IQR(diff(y))/1.908

  #
  # first check arguments and initialize
  #
  # test dimension of data (vector of 3D) and define dimension related stuff
  if(is.null(res)) {
      return(warning("Please specify keep=''all'' when calling fmri.lm"))
  }
  d <- 3
  dy <- dim(y)
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1*n2*n3
  if (length(dy)==d) {
    dim(y) <- dy <- c(dy,1)
  } else if (length(dy)!=d+1) {
    stop("y has to be 3 or 4 dimensional")
  }
  dv <- dim(y)[d+1]
  # define vwghts
  if (length(vwghts)>dv) vwghts <- vwghts[1:dv]
  dv0 <- sum(vwghts>0)

  # MAE
  mae <- NULL

  # set the code for the kernel (used in lkern) and set lambda
  lkern <- switch(lkern,
                  Triangle=2,
                  Plateau=1,
                  Gaussian=3,
                  1)
  skern <- switch(skern,
                  Triangle=2,
                  Plateau=1,
                  Exponential=3,
                  2)

  # define qlambda, lambda
  if (is.null(qlambda)) qlambda <- switch(skern,.9945,.9988,.9995)
  if (qlambda<.9) warning("Inappropriate value of qlambda")
  if (qlambda<1) {
    lambda <- qchisq(qlambda,dv0) 
  } else {
    lambda <- 1e50
  }
  if(skern%in%c(1,2)) {
    # to have similar preformance compared to skern="Exp"
    lambda <- 4/3*lambda
    if(skern==1) if(is.null(spmin)) spmin <- .3 
    spmax <- 1
  } else {
      spmin <- 0
      spmax <- 4
  }
  # set hinit and hincr if not provided
  if (is.null(hinit)||hinit<1) hinit <- 1
  
  # define hmax
  if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points

  # re-define bandwidth for Gaussian lkern!!!!
  if (lkern==3) {
    # assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- fwhm2bw(hmax)*4
    hinit <- min(hinit,hmax)
  }
  if (qlambda == 1) hinit <- hmax

  # define hincr
  if (is.null(hincr)||hincr<=1) hincr <- 1.25
  hincr <- hincr^(1/d)

  if (length(dim(sigma2))!=3) stop("insufficient variance information for full analysis") 
  if(is.null(mask)) mask <- array(TRUE,dy[1:3])
  mask[sigma2>=1e16] <- FALSE
#  in these points sigma2 probably contains NA's
  # deal with homoskedastic Gaussian case by extending sigma2
  if (length(sigma2)!=n) stop("sigma2 does not have length 1 or same length as y")
  dim(sigma2) <- dy[1:3]
  sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 

  # demo mode should deliver graphics!
  if (demo && !graph) graph <- TRUE
  
  # Initialize  list for bi and theta
  if (is.null(wghts)) wghts <- c(1,1,1)
  hinit <- hinit/wghts[1]
  hmax <- hmax/wghts[1]
  wghts <- (wghts[2:3]/wghts[1])
  tobj <- list(bi= sigma2)
  theta <- y
  maxvol <- getvofh(hmax,lkern,wghts)
  kstar <- as.integer(log(maxvol)/log(1.25))
  steps <- kstar+1
  
  if(qlambda < 1) k <- 1 else k <- kstar
  hakt <- hinit
  hakt0 <- hinit
  lambda0 <- lambda 
  if (hinit>1) lambda0 <- 1e50 # that removes the stochstic term for the first step
  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  if(testprop) {
    if(is.null(u)) u <- 0
  } 
  propagation <- NULL
  # run single steps to display intermediate results
  residuals <- readBin(res,"integer",prod(ddim),2)*resscale
  dim(residuals) <- c(ddim[4],ddim[1:3])
  vadjust <- apply(residuals,2:4,var)*sigma2
  while (k<=kstar) {
      hakt0 <- gethani(1,3,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,3,lkern,1.25^k,wghts,1e-4)
      hakt.oscale <- if(lkern==3) bw2fwhm(hakt/4) else hakt
      cat("step",k,"bandwidth",signif(hakt.oscale,3)," ")
    dlw <- (2*trunc(hakt/c(1,wghts))+1)[1:d]
    hakt0 <- hakt
    theta0 <- theta
    bi0 <- tobj$bi/(1+log(tobj$bi/sigma2))
    #
    #   need these values to compute variances after the last iteration
    #
    tobj <- .Fortran("chawsv",as.double(y),
                     as.double(residuals),
                     as.double(sigma2),
                     as.logical(!mask),
                     as.logical(weighted),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.integer(ddim[4]),
                     as.integer(dv),
                     as.integer(dv0),
                     hakt=as.double(hakt),
                     as.double(lambda0),
                     as.double(theta0),
                     bi=as.double(bi0),
                     resnew=double(prod(ddim)),
                     thnew=double(n1*n2*n3*dv),
                     as.integer(lkern),
                     as.integer(skern),
                     as.double(spmin),
                     as.double(spmax),
                     double(prod(dlw)),
                     as.double(wghts),
                     as.double(vwghts),
                     double(dv),#swjy
                     double(dv0),#thi
                     double(ddim[4]),#resi
                     PACKAGE="fmri",DUP=TRUE)[c("bi","thnew","hakt","resnew")]
    gc()
    theta <- array(tobj$thnew,dy) 
    dim(tobj$bi) <- ddim[1:3]
    tobj$bi <- tobj$bi*vadjust
#  correcting for missing term in variance estimate from residuals, e.g. effects of design matrix
#  now contains 1/var(thetahat)
    if(testprop) {
      pobj <- .Fortran("chawsv",as.double(y),
                        as.double(residuals),
                        as.double(sigma2),
                        as.logical(!mask),
                        as.logical(weighted),
                        as.integer(n1),
                        as.integer(n2),
                        as.integer(n3),
                        as.integer(ddim[4]),
                        as.integer(dv),
                        as.integer(dv0),
                        hakt=as.double(hakt),
                        as.double(1e50),
                        as.double(theta0),
                        bi=as.double(bi0),
                        resnew=double(prod(ddim)),
                        thnew=double(n1*n2*n3*dv),
                        as.integer(lkern),
                        as.integer(skern),
                        as.double(spmin),
                        as.double(spmax),
                        double(prod(dlw)),
                        as.double(wghts),
                        as.double(vwghts),
                        double(dv),#swjy
                        double(dv0),#thi
                        double(ddim[4]),#resi
		        PACKAGE="fmri",DUP=TRUE)[c("bi","thnew","resnew")]
      ptheta <- array(pobj$thnew,dy) 
      rm(pobj) 
      gc()
      propagation <- c(propagation,sum(abs(theta-ptheta))/sum(abs(ptheta-u)))
      cat("Propagation with alpha=",max(propagation),"\n")
      cat("alpha values:","\n")
      print(rbind(hakt.oscale,signif(propagation[-1],3)))
    }
    if (graph) {
      par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
      image(y[,,n3%/%2+1,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
      title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
      image(theta[,,n3%/%2+1,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
      title(paste("Reconstruction  h=",signif(hakt.oscale,3)," min=",signif(min(theta),3)," max=",signif(max(theta),3)))
      image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
      title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
    }
    if (!is.null(u)) {
      cat("bandwidth: ",signif(hakt.oscale,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
          signif(mean((theta-u)^2),3),"   MAE: ",signif(mean(abs(theta-u)),3)," mean(bi)=",signif(mean(tobj$bi),3),"\n")
      mae <- c(mae,signif(mean(abs(theta-u)),3))
    } else if (max(total) >0) {
      cat(signif(total[k],2)*100,"%                 \r",sep="")
     }
    if (demo) readline("Press return")
    k <- k+1
#  adjust lambda for the high intrinsic correlation between  neighboring estimates 
#    quot <- .Fortran("ni2var",as.double(hakt),
#                              as.integer(lkern),
#                              as.double(wghts),
#                              quot=double(1))$quot
#    cat("ni/var",quot,"\n")
                              
    lambda0 <- if(k<=length(ladjust)) ladjust[k]*lambda else ladjust*lambda#/quot
    gc()
  }
  cat("fmri.smooth: estimate correlations","\n")
  lags <- c(5,5,3)
  scorr <- .Fortran("imcorr",as.double(tobj$resnew),
                     as.logical(mask),
                     as.integer(n3),
                     as.integer(n2),
                     as.integer(n1),
                     as.integer(ddim[4]),
                     scorr=double(prod(lags)),
                     as.integer(lags[1]),
                     as.integer(lags[2]),
                     as.integer(lags[3]),
                     PACKAGE="fmri",DUP=TRUE)$scorr
  dim(scorr) <- lags                     
  vartheta <- 1/tobj$bi
  vred <- vartheta*sigma2
  z <- list(theta=theta,ni=tobj$bi,var=vartheta,vred=vred,vred0=median(vred[mask]),y=y,
            hmax=tobj$hakt*switch(lkern,1,1,bw2fwhm(1/4)),mae=mae,
            alpha=propagation,scorr=scorr,res=tobj$resnew,mask=mask)
  #   vred accounts for variance reduction with respect to uncorrelated (\check{sigma}^2) data
  class(z) <- "aws.gaussian"
  invisible(z)
}


smooth3D <- function(y,lkern="Gaussian",weighted=FALSE,sigma2=NULL,mask=NULL,hmax=NULL,
                     wghts=NULL) {
  #
  # first check arguments and initialize
  #
  # test dimension of data (vector of 3D) and define dimension related stuff
  d <- 3
  dy <- dim(y)
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1*n2*n3
  if (length(dy)==d) {
    dim(y) <- dy <- c(dy,1)
  } else if (length(dy)!=d+1) {
    stop("y has to be 3 or 4 dimensional")
  }
  dv <- dim(y)[d+1]
    if(is.null(sigma2)) {
       weighted <- FALSE
    } else {
      if(length(sigma2)!=n) weighted <- FALSE
      sigma2 <- 1/sigma2
    }
  if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points

  # re-define bandwidth for Gaussian lkern!!!!
  lkern <- switch(lkern,
                  Triangle=2,
                  Plateau=1,
                  Gaussian=3,
                  1)
  if (lkern==3) {
    # assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- fwhm2bw(hmax)*4
  }
  if (is.null(wghts)) wghts <- c(1,1,1)
  if(is.null(mask)) mask <- array(TRUE,dy[1:3])
  hmax <- hmax/wghts[1]
  wghts <- (wghts[2:3]/wghts[1])
    dlw <- (2*trunc(hmax/c(1,wghts))+1)[1:d]
    ysmooth <- .Fortran("smooth3d",
                     as.double(y),
                     as.double(sigma2),
                     as.logical(!mask),
                     as.logical(weighted),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.integer(dv),
                     hakt=as.double(hmax),
                     thnew=double(n1*n2*n3*dv),
                     as.integer(lkern),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(dv),#swjy
                     PACKAGE="fmri",DUP=TRUE)$thnew
array(ysmooth,dy)
}
