###########################################################################
#
#   nonadaptive 1D -- 3D smoothing on a grid (Gaussian kernel)
#
###########################################################################
kernsm<-function (y, h = 1, kern="Gaussian", m=0, nsector=1, sector=1, symmetric=FALSE)
{
#
#  nonadaptive kernel smoothing using FFT
#
    args <- match.call()
    expand.x.par <- function(x,h){
       dx <- dim(x)
       if(is.null(dx)) dx <- length(x)
       d <- length(dx)
       if(length(h)<d) h <- pmax(.01,rep(h[1],d))
#      zero-padding
       dx1 <- nextn(dx+2*h)
       ddx <- (dx1-dx)%/%2
       ilow <- ddx+1
       iup <- ddx+dx
       list(dx1=dx1,d=d,ilow=ilow,iup=iup,h=h)
    }
    expand.x <- function(x,xp){
       xx <- array(0,xp$dx1)
       if(xp$d==1) {
          xx[xp$ilow:xp$iup] <- x 
       } else if(xp$d==2) {
          xx[xp$ilow[1]:xp$iup[1],xp$ilow[2]:xp$iup[2]] <- x 
       } else {
          xx[xp$ilow[1]:xp$iup[1],xp$ilow[2]:xp$iup[2],xp$ilow[3]:xp$iup[3]] <- x 
       }
       xx
    }
    grid <- function(d) {
       d0 <- d%/%2+1
       gd <- seq(0,1,length=d0)
       if (2*d0==d+1) gd <- c(gd,-gd[d0:2]) else gd <- c(gd,-gd[(d0-1):2])
       gd
    }
    gridind <- function(d,side=1) {
       d0 <- d%/%2+1
       if(side==1) 1:d0  else (d0+1):d 
    }
    lkern1 <- function(x,h,ikern,m){
       nx <- length(x)
       .Fortran("lkern1",
                as.double(x),
                as.integer(nx),
                as.double(h),
                as.integer(ikern),
                as.integer(m),
                khofx=double(nx),
                DUPL=TRUE,
                PACKAGE="aws")$khofx
    }
    lkern <- function(xp,kind="Gaussian",m=0,nsector=1,sector=1,symmetric=FALSE){
#
#   generate 1D, 2D or 3D product kernel weights appropriate for fft
#   m defines (partial) derivatives up to order 2
#   
#
       xh <- array(0,c(xp$d,xp$dx1))
       dx0 <- xp$dx1%/%2 + 1
       if(length(m)<xp$d) m <- rep(m,xp$d)
       if(xp$d==1) {
          xh[1,] <- grid(xp$dx1)
       } else if(xp$d==2) {
          xh[1,,] <- grid(xp$dx1[1])
          for(i in 1:xp$dx1[1]) xh[2,i,] <- grid(xp$dx1[2])
       } else if(xp$d==2) {
          xh[1,,,] <- grid(xp$dx1[1])
          for(i in 1:xp$dx1[1]) xh[2,i,,] <- grid(xp$dx1[2])
          for(i in 1:xp$dx1[1]) for(j in 1:xp$dx1[2]) xh[2,i,j,] <- grid(xp$dx1[3])
       }
       ikern <- switch(kind,"Gaussian"=1,
                            "Uniform"=2,
                            "Triangle"=3,
                            "Epanechnikov"=4,
                            "Biweight"=5,
                            "Triweight"=6,
                            1)
# use  equivalent kernels for local polynomial smoothing with p=m+1                
       if(any(m>2)) {
          m <- pmin(m,2)
          warning("only second order derivatives implemented, setting m==2")
       }
       kwghts <- switch(xp$d,lkern1(grid(xp$dx1[1]),2*xp$h[1]/xp$dx1[1],ikern,m[1]),
                 outer(lkern1(grid(xp$dx1[1]),2*xp$h[1]/xp$dx1[1],ikern,m[1]),
                       lkern1(grid(xp$dx1[2]),2*xp$h[2]/xp$dx1[2],ikern,m[2]),"*"),
                 outer(outer(lkern1(grid(xp$dx1[1]),2*xp$h[1]/xp$dx1[1],ikern,m[1]),
                       lkern1(grid(xp$dx1[2]),2*xp$h[2]/xp$dx1[2],ikern,m[2]),"*"),
                       lkern1(grid(xp$dx1[3]),2*xp$h[3]/xp$dx1[3],ikern,m[3]),"*"))
       if(nsector>1){
          if(xp$d==1){
             if(sector%in%c(1,2)){
                d0 <- xp$dx1%/%2+1
                if(sector==2){
                   kwghts[(d0+1):xp$dx1] <- 0
                } else {
                   kwghts[1:(d0-1)] <- 0 
                }
             } else {
                warning("sector needs to be 1 (left sided) or 2 (right sided)") 
             }
          }
          if(xp$d==2){
             sector <- .Fortran("sector",
                                as.double(grid(xp$dx1[1])),
                                as.integer(xp$dx1[1]),
                                as.double(grid(xp$dx1[2])),
                                as.integer(xp$dx1[2]),
                                as.integer(nsector),
                                as.integer(sector),
                                as.logical(symmetric),
                                insector=double(xp$dx1[1]*xp$dx1[2]),
                                DUPL=TRUE,
                                PACKAGE="aws")$insector                               
             kwghts <- kwghts*array(sector,xp$dx1)
          }
          if(xp$d==3){
             sect3D <- switch(sector,c(0,0,0),
                                     c(0,0,1),
                                     c(0,1,0),
                                     c(0,1,1),
                                     c(1,0,0),
                                     c(1,0,1),
                                     c(1,1,0),
                                     c(1,1,1)) 
             sector <- array(0,xp$dx1)
             sector[gridind(xp$dx1[1]),
                    gridind(xp$dx1[2]),
                    gridind(xp$dx1[3])] <- 1
             if(symmetric) sector[-gridind(xp$dx1[1]),
                                  -gridind(xp$dx1[2]),
                                  -gridind(xp$dx1[3])] <- 1

             kwghts <- kwghts*array(sector,xp$dx1)
          }
          kwghts <- kwghts/sum(kwghts)
       }
       kwghts
    }
    ypar <- expand.x.par(y,h)
    yext <- expand.x(y,ypar)
    kwghts <- lkern(ypar,kern,m,nsector,sector,symmetric)
    yhat <- Re(fft(fft(yext) * fft(kwghts),inverse=TRUE))/prod(ypar$dx1)
    ilow <- ypar$ilow
    iup <- ypar$iup
    yhat <- switch(ypar$d,yhat[ilow:iup],
                          yhat[ilow[1]:iup[1],ilow[2]:iup[2]],
                          yhat[ilow[1]:iup[1],ilow[2]:iup[2],ilow[3]:iup[3]])
    kernsmobj(y,h=h,kern=kern,m=m,nsector=nsector,sector=sector,symmetric=symmetric,
              yhat=yhat,vred=1/sum(kwghts^2),call=args)
    }
ICIsmooth <- function(y, hmax, hinc=1.45, thresh=NULL, kern="Gaussian", m=0, sigma=NULL, nsector=1, 
             sector=1, symmetric=FALSE, presmooth = FALSE){
   args <- match.call()
   if(any(m>0)&nsector>1){
      nsector <- 1
      sector <- 1
      warning("no sectorial weights for estimates of derivatives")
   }
   beta <- .05
   n <- length(y)
   dy <- dim(y)
   if(is.null(dy)) dy <- n
   d <- length(dy)
   if(d==1&nsector>2-symmetric){
      warning(paste("ICIsmooth: only nsector%in%c(1,",(2-symmetric),") implemented in 1D\n 
                    using nsector=",2-symmetric,sep=""))
      nsector <- (2-symmetric)
   } 
   if(d>2&nsector>1&nsector!=4*(2-symmetric)){
      warning(paste("ICIsmooth: nsector%in%c(1,",4*(2-symmetric),") implemented in 3D\n 
            using nsector=",4*(2-symmetric),sep=""))
      nsector <- 4*(2-symmetric)
   }
   if(sector<1||sector>nsector){
      warning("ICIsmooth: illegal value of sector, setting sector=1")
   }
   if(length(m)<d) m <- rep(m,d)
   if(is.null(thresh)) {
      thresh <- sqrt((d+2*sum(m))/2)+qnorm(1-beta/2)
      cat("using thresh=",thresh,"\n")
      }
   if(is.null(sigma)) sigma <- median(abs(y[-1]-y[-n]))/.9538
   if(all(m==0)){
   Low <- as.vector(y - thresh * sigma)
   Up  <- as.vector(y + thresh * sigma)
   hakt <- if(kern=="Gaussian") .3 else 1
   } else {
   Low <- rep(-Inf,n)
   Up  <- rep(Inf,n)
   hakt <- if(kern=="Gaussian") .6 else hinc*(max(m)+1)
   }
   hbest <- rep(hakt,n)
   fixed <- rep(FALSE,n)
   yhat <- as.vector(y)
   vhat <- rep(sigma^2,n) 
   while(hakt < hmax){
      z <- kernsm(y, hakt, kern, m, nsector, sector, symmetric)
      ind0 <- (1:n)[!fixed]
      Low[ind0] <- pmax(Low,z@yhat-thresh * sigma/sqrt(z@vred))[ind0]
      Up[ind0] <- pmin(Up,z@yhat+thresh * sigma/sqrt(z@vred))[ind0]
      ind <- ind0[Low[ind0]<=Up[ind0]]
      hbest[ind] <- hakt
      yhat[ind] <- z@yhat[ind]
      vhat[ind] <- sigma^2/z@vred
      fixed[-ind] <- TRUE
      hakt <- hakt*hinc
      if(sum(fixed)==n) break
   } 
   if(presmooth){
      hbest <- switch(d,.Fortran("median1d",
                                 as.double(hbest),
                                 as.integer(n),
                                 hbest=double(n),
                                 DUPL=TRUE,
                                 PACKAGE="aws")$hbest,
                        .Fortran("median2d",
                                 as.double(hbest),
                                 as.integer(dy[1]),
                                 as.integer(dy[2]),
                                 hbest=double(n),
                                 DUPL=TRUE,
                                 PACKAGE="aws")$hbest,
                        .Fortran("median3d",
                                 as.double(hbest),
                                 as.integer(dy[1]),
                                 as.integer(dy[2]),
                                 as.integer(dy[3]),
                                 hbest=double(n),
                                 DUPL=TRUE,
                                 PACKAGE="aws")$hbest)
      hakt <- if(all(m==0)) { if(kern=="Gaussian") .3 else 1 }
      else { if(kern=="Gaussian") .6 else hinc*(max(m)+1) }
      while(hakt < hmax){
         ind <- (1:n)[abs(hakt-hbest)<1e-3]
         if(length(ind)>0){
            z <- kernsm(y, hakt, kern, m, nsector, sector, symmetric)
            yhat[ind] <- z@yhat[ind]
            vhat[ind] <- sigma^2/z@vred
         }
         hakt <- hakt*hinc
      }                                  
   }
   ICIsmoothobj(y,h=hmax,hinc=hinc,thresh=thresh,kern=kern,m=m,nsector=nsector,
                sector=sector,symmetric=symmetric,yhat=array(yhat,dy),
                vhat=vhat,hbest=hbest,sigma=sigma,call=args)
}
ICIcombined <- function(y, hmax, hinc=1.45, thresh=NULL, kern="Gaussian", m=0, 
                 sigma = NULL, nsector=1, symmetric=FALSE, presmooth=FALSE){
   args <- match.call()
   if(any(m>0)&nsector>1){
      nsector <- 1
      warning("no sectorial weights for estimates of derivatives")
   }
   n <- length(y)
   dy <- dim(y)
   if(is.null(dy)) dy <- n
   d <- length(dy)
   if(d==1&nsector>2-symmetric){
      warning(paste("ICIcombined: only nsector%in%c(1,",(2-symmetric),") implemented in 1D\n 
                    using nsector=",2-symmetric,sep=""))
      nsector <- (2-symmetric)
   } 
   if(d>2&nsector>1&nsector!=4*(2-symmetric)){
      warning(paste("ICIcombined: nsector%in%c(1,",4*(2-symmetric),") implemented in 3D\n 
            using nsector=",4*(2-symmetric),sep=""))
      nsector <- 4*(2-symmetric)
   }
   if(is.null(sigma)) sigma <- median(abs(y[-1]-y[-n]))/.9538
   if(nsector==1){
      return(ICIsmooth(y, hmax, hinc, thresh, kern, m, 1, 1, symmetric, presmooth))
   } else {
      yhatc <- array(0,c(nsector,prod(dy)))
      vhatc <- array(0,c(nsector,prod(dy)))
      hbest <- numeric(prod(dy))
      for(i in 1:nsector){
         z <- ICIsmooth(y, hmax, hinc, thresh, kern, m, sigma, nsector, i, symmetric, presmooth)
         yhatc[i,] <- z@yhat
         vhatc[i,] <- z@vhat
         hbest <- hbest+z@hbest
      }
      hbest <- hbest/nsector
      vhatinv <- 1/vhatc[1,]
      for(i in 2:nsector) vhatinv <- vhatinv+1/vhatc[i,]
      vhat <- 1/vhatinv
      yhat <- vhat/vhatc[1,]*yhatc[1,]
      for(i in 2:nsector) yhat <- yhat+vhat/vhatc[i,]*yhatc[i,]
   }
   ICIsmoothobj(y,h=hmax,hinc=hinc,thresh=thresh,kern=kern,m=m,nsector=nsector,
                sector=0,symmetric=symmetric,yhat=array(yhat,dy),vhat=vhat,
                hbest=hbest,sigma=sigma,call=args)
}
riskyhat <- function(y,u){
if(length(u)==1) {
u <- rep(u,length(y))
dim(u) <- dim(y)
}
if(is.null(dim(y))||length(dim(y))==1) {
   if(length(y)!=length(u)) stop("length(y)!=length(u)") else u <- as.vector(u)
} else {
   if(length(dim(y))!=length(dim(u))||any(dim(y)!=dim(u))) stop("dim(y)!=dim(u)")
}
MSE <- mean((y-u)^2)
RMSE <- sqrt(MSE)
SNR <- 10*log(mean(u^2)/MSE,10)
PSNR <- 10*log(max(u^2)/MSE,10)
MAE <- mean(abs(y-u))
MaxAE <- max(abs(y-u))
ymean <- mean(y)
umean <- mean(u)
sigy <- sd(as.vector(y))
sigu <- sd(as.vector(u))
# now the Universal image quality index of Wang and Bovik (2002)
UIQI <- cor(as.vector(y),as.vector(u))*2*umean*ymean/(umean^2+ymean^2)*
        2*sigy*sigu/(sigy^2+sigu^2)
cat("RMSE=",signif(RMSE,4)," SNR=",signif(SNR,3)," MAE=",signif(MAE,4),
    " MaxAE=",signif(MaxAE,3)," UIQI=",signif(UIQI,4),"\n")
invisible(list(RMSE=RMSE,SNR=SNR,PSNR=PSNR,MAE=MAE,MaxAE=MaxAE,UIQI=UIQI))
}


