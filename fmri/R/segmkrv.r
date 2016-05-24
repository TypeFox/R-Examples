#########################################################################################################################
#
#    R - function  segm3Dkrv  for simulating critical values in segm3D
#
#    emaphazises on the propagation-separation approach 
#
#    Copyright (C) 2010-12 Weierstrass-Institut fuer
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

segm3Dkrv <- function(dy,df,hmax=NULL,ladjust=1,beta=0,graph=FALSE,h0=c(0,0,0)) {
#
#
#  Auxilary functions
   IQRdiff <- function(y) IQR(diff(y))/1.908
#
# first check arguments and initialize
#
   args <- match.call()
   nt <- df+1
   if (length(dy)!=3) {
      stop("dy has to be of length 3")
   }
   d <- 3
   n1 <- dy[1]
   n2 <- dy[2]
   n3 <- dy[3]
   n <- n1*n2*n3
   res <- array(rnorm(prod(dy)*nt),c(nt,dy))
   if(any(h0>0)) {
#      require(aws)
      warning("for simulating critical values we need package aws")
#      for(i in 1:nt) res[i,,,] <- kernsm(res[i,,,],h0)@yhat
   }
# test dimension of data (vector of 3D) and define dimension related stuff
   ddim <- dim(res)
   y <- .Fortran("mean3D",
                 as.double(res),
                 as.integer(n1),
                 as.integer(n2),
                 as.integer(n3),
                 as.integer(nt),
                 y=double(prod(dy)),
                 PACKAGE="fmri",DUP=TRUE)$y
   dim(y) <- dy
   if (length(dy)==d+1) {
      dim(y) <- dy[1:3]
   } else if (length(dy)!=d) {
      stop("y has to be 3 dimensional")
   }
# set the code for the kernel (used in lkern) and set lambda
   lkern <- 1
   skern <- 1
# define lambda
   lambda <- ladjust*(exp(2.6-3.17*log(df)+8.4*log(log(df)))+16) # corresponding to p_0 ~ 1e-6
   hinit <- 1
  
# define hmax
   if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points

# re-define bandwidth for Gaussian lkern!!!!
   if (lkern==3) {
# assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
      hmax <- fwhm2bw(hmax)*4
      hinit <- min(hinit,hmax)
   }
   if(is.null(h0)) h0 <- rep(0,3)
# estimate variance in the gaussian case if necessary  
# deal with homoskedastic Gaussian case by extending sigma2
  mask <- array(TRUE,dy[1:3])
  res <- .Fortran("sweepm",res=as.double(res),
                           as.logical(mask),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(n3),
                           as.integer(nt),
                           PACKAGE="fmri",DUP=TRUE)$res                         
  cat("\nfmri.smooth: first variance estimate","\n")
  vartheta0 <- .Fortran("ivar",as.double(res),
                           as.double(1),
                           as.logical(rep(TRUE,prod(dy))),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(n3),
                           as.integer(nt),
                           var = double(n1*n2*n3),
                           PACKAGE="fmri",DUP=TRUE)$var
   sigma2 <- vartheta0/df #  thats the variance of  y  ... !!!! assuming zero mean
   sigma2 <- 1/sigma2 # need the inverse for easier computations
   dim(sigma2) <- dy
# Initialize  list for bi and theta
   wghts <- c(1,1,1)
   hinit <- hinit/wghts[1]
   hmax <- hmax/wghts[1]
   wghts <- (wghts[2:3]/wghts[1])
   tobj <- list(bi= rep(1,n))
   theta <- y
   segm <- array(0,dy)
   varest <- 1/sigma2
   maxvol <- getvofh(hmax,lkern,wghts)
   fov <- prod(ddim[1:3])
   kstar <- as.integer(log(maxvol)/log(1.25))
   steps <- kstar+1
   cat("FOV",fov,"ladjust",ladjust,"lambda",lambda,"\n")
   k <- 1 
   hakt <- hinit
   hakt0 <- hinit
   lambda0 <- lambda
   maxvalue <- matrix(0,2,kstar)
   mse <- numeric(kstar)
   mae <- numeric(kstar)
   if (hinit>1) lambda0 <- 1e50 # that removes the stochstic term for the first step
   scorr <- numeric(3)
   if(h0[1]>0) scorr[1] <-  get.corr.gauss(h0[1],2)
   if(h0[2]>0) scorr[2] <-  get.corr.gauss(h0[2],2)
   if(h0[3]>0) scorr[3] <-  get.corr.gauss(h0[3],2)
   total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
# run single steps to display intermediate results
   while (k<=kstar) {
      hakt0 <- gethani(1,10,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,10,lkern,1.25^k,wghts,1e-4)
      hakt.oscale <- if(lkern==3) bw2fwhm(hakt/4) else hakt
      cat("step",k,"bandwidth",signif(hakt.oscale,3)," ")
      dlw <- (2*trunc(hakt/c(1,wghts))+1)[1:d]
      hakt0 <- hakt
      theta0 <- theta
      bi0 <- tobj$bi
#
#   need these values to compute variances after the last iteration
#
      tobj <- .Fortran("segm3dkb",
                       as.double(y),
                       as.double(res),
                       as.double(sigma2),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.integer(nt),
                       as.double(df),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(theta0),
                       bi=as.double(bi0),
                       thnew=double(n1*n2*n3),
                       as.integer(lkern),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nt),#swres
                       as.double(fov),
                       varest=as.double(varest),
                       maxvalue=double(1),
                       minvalue=double(1),
                       PACKAGE="fmri",DUP=TRUE)[c("bi","thnew","hakt","varest","maxvalue","minvalue")]
      gc()
      theta <- array(tobj$thnew,dy) 
      varest <- array(tobj$varest,dy)
      dim(tobj$bi) <- dy
      maxvalue[1,k] <- tobj$maxvalue
      maxvalue[2,k] <- -tobj$minvalue
      mae[k] <- mean(abs(theta))
      mse[k] <- mean(theta^2)
      if (graph) {
         par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
         image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
         image(theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Reconstruction  h=",signif(hakt.oscale,3)," min=",signif(min(theta),3),"   max=",signif(max(theta),3)))
         image(segm[,,n3%/%2+1]>0,col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Segmentation  h=",signif(hakt.oscale,3)," detected=",sum(segm>0)))
         image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
      }
      if (max(total) >0) {
         cat(signif(total[k],2)*100,"%                 \r",sep="")
      }
      k <- k+1
#  adjust lambda for the high intrinsic correlation between  neighboring estimates 
      lambda0 <- lambda
      gc()
   }

  z <- list(mae=mae,mse=mse,maxvalue=maxvalue)
  invisible(z)
}
