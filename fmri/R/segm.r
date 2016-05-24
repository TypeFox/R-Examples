###################################################################################################
#
#    R - function  segm3D  for Adaptive Weights Segmentation in 3D SPM's
#    
#    exact value for Variance reduction from smoothed residuals
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

segm3D <- function(y,weighted=TRUE,
                   sigma2=NULL,mask=NULL,hinit=NULL,hmax=NULL,
                   ladjust=1,graph=FALSE,wghts=NULL,
                   df=100,h0=c(0,0,0),res=NULL, resscale=NULL, 
                   ddim=NULL,delta=0,alpha=.05,restricted=FALSE) {
#
#
#  Auxilary functions
   IQRdiff <- function(y) IQR(diff(y))/1.908
   getkrval <- function(df,ladj,fov,k,alpha){
      dfinv <- 1/df
      lfov <- log(fov)
      idffov <- dfinv*log(fov)
      lalpha <- log(alpha)
      lk <- log(k)
      ladj <- ladj-1
      kv <- 1.0105241 + 3.2036997*dfinv + 0.5442412*idffov -
            0.0002686*lalpha -0.0172885*ladj +0.0009396*lk -1.7976621*dfinv*log(df) -
            0.7856579*dfinv*lalpha+0.0765149*idffov*ladj
      kv
      }
#
# first check arguments and initialize
#
   args <- match.call()
   if(is.null(ddim)) ddim <- dim(y)
# test dimension of data (vector of 3D) and define dimension related stuff
   d <- 3
   dy <- dim(y)
   n1 <- dy[1]
   n2 <- dy[2]
   n3 <- dy[3]
   n <- n1*n2*n3
   nt <- ddim[4]
   if (length(dy)==d+1) {
      dim(y) <- dy[1:3]
   } else if (length(dy)!=d) {
      stop("y has to be 3 dimensional")
   }
# set the code for the kernel (used in lkern) and set lambda
   lkern <- 1
   skern <- 1
# define lambda
   lambda <- ladjust*(exp(2.6-3.17*log(df)+8.4*log(log(df)))+16)
# to have similar preformance compared to skern="Exp"
   if (is.null(hinit)||hinit<1) hinit <- 1
# define hmax
   if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points
# estimate variance in the gaussian case if necessary  
# deal with homoskedastic Gaussian case by extending sigma2
   if (length(sigma2)==1) sigma2<-array(sigma2,dy[1:3]) 
   if (length(sigma2)!=n) stop("sigma2 does not have length 1 or same length as y")
   dim(sigma2) <- dy[1:3]
   if(is.null(mask)) mask <- array(TRUE,dy[1:3])
   mask[sigma2>=1e16] <- FALSE
#  in these points sigma2 probably contains NA's
   sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
# deal with homoskedastic Gaussian case by extending sigma2
   residuals <- readBin(res,"integer",prod(ddim),2)
  cat("\nfmri.smooth: first variance estimate","\n")
  varest0 <- .Fortran("ivar",as.double(residuals),
                           as.double(resscale),
                           as.logical(mask),
                           as.integer(n1),
                           as.integer(n2),
                           as.integer(n3),
                           as.integer(nt),
                           var = double(n1*n2*n3),
                           PACKAGE="fmri",DUP=TRUE)$var
   vq <- varest0*sigma2
#   plot(density(vq),main="Density of vq")
   if (is.null(wghts)) wghts <- c(1,1,1)
   hinit <- hinit/wghts[1]
   hmax <- hmax/wghts[1]
   wghts <- wghts[2:3]/wghts[1]
   tobj <- list(bi= sigma2)
   theta <- y
   segm <- array(0,dy[1:3])
   varest <- varest0
   maxvol <- getvofh(hmax,lkern,wghts)
   fov <- sum(mask)
   kstar <- as.integer(log(maxvol)/log(1.25))
   steps <- kstar+1
   k <- 1 
   hakt <- hinit
   hakt0 <- hinit
   lambda0 <- lambda
   if (hinit>1) lambda0 <- 1e50 # that removes the stochstic term for the first step
   scorr <- numeric(3)
   total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
#   for(i in 10:kstar) thresh <- max(thresh,tadj*getkrval(df,ladjust,fov,i,alpha))
   thresh <- getkrval(df,ladjust,fov,kstar,alpha)
#  just to ensure monotonicity of thresh with kmax, there exist a few parameter configurations
#  where the approximation formula does not ensure monotonicity
   cat("FOV",fov,"delta",delta,"thresh",thresh,"ladjust",ladjust,"lambda",lambda,"df",df,"\n")
   residuals <- residuals*resscale
   kv <- array(0,dy[1:3])
#
#   need these values to compute variances 
#
   while (k<=kstar) {
      hakt0 <- gethani(1,10,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,10,lkern,1.25^k,wghts,1e-4)
      cat("step",k,"bandwidth",signif(hakt,3)," ")
      dlw <- (2*trunc(hakt/c(1,wghts))+1)[1:d]
      hakt0 <- hakt
      bi0 <- tobj$bi
      tobj <- .Fortran("segm3d",
                       as.double(y),
                       as.double(residuals),
                       as.double(sigma2),
                       as.logical(!mask),
                       as.logical(weighted),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.integer(nt),
                       as.double(df),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(theta),
                       bi=as.double(bi0),
                       thnew=double(n1*n2*n3),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nt),#swres
                       double(n1*n2*n3),#pvalues
                       segm=as.integer(segm),
                       as.double(delta),
                       as.double(thresh),
                       as.double(fov),
                       as.double(vq),
                       as.double(varest0),
                       varest=as.double(varest),
                       as.logical(restricted),
                       PACKAGE="fmri",DUP=TRUE)[c("bi","thnew","hakt","segm","varest")]
      gc()
      theta <- array(tobj$thnew,dy[1:3]) 
      segm <- array(tobj$segm,dy[1:3])
      varest <- array(tobj$varest,dy[1:3])
      dim(tobj$bi) <- dy[1:3]
      if (graph) {
         par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
         image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
         image(theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(theta),3),"   max=",signif(max(theta),3)))
         image(segm[,,n3%/%2+1]>0,col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Segmentation  h=",signif(hakt,3)," detected=",sum(segm>0)))
         image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
      }
      if (max(total) >0) {
         cat(signif(total[k],2)*100,"%                 \r",sep="")
      }
      k <- k+1
      lambda0 <- lambda
      gc()
   }
  z <- list(theta=theta,ni=tobj$bi,var=varest,y=y,segm=segm,
            hmax=tobj$hakt*switch(lkern,1,1,bw2fwhm(1/4)),
            scorr=scorr,mask=mask)
  class(z) <- "aws.gaussian"
  invisible(z)
}

