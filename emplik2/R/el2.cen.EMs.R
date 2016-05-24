el2.cen.EMs<-function(x,dx,y,dy,fun=function(x,y) {x>=y}, mean=0.5,
  maxit=25){

#x,y pairs can be any combination of uncensored, left-cens, right-cens
#Data can be discrete or continuous
#Note that x>y is not the same as x>=y for discrete data

#Store data x,y as vectors
    xvec <- as.vector(x)
    yvec <- as.vector(y)

#Store length of xvec,yvec   
    nx <- length(xvec)
    ny <- length(yvec)

#Check that there are at least 2 data in x,y
    if (nx <= 1) 
        stop("need more observations in x")
    if (ny <= 1) 
        stop("need more observations in y")

#Check that status and observations have same length in x,y
    if (length(dx) != nx) 
        stop("length of x and dx must agree")
    if (length(dy) != ny) 
        stop("length of y and dy must agree")

#Check that status are only 0,1,2
    if (any((dx != 0) & (dx != 1) & (dx != 2)))  
        stop("dx must be 0(right-censored) or 1(uncensored) or 2(left-
         censored")
    if (any((dy != 0) & (dy != 1) & (dy != 2)))  
        stop("dy must be 0(right-censored) or 1(uncensored) or 2(left-
         censored") 
   
#Check that xvec,yvec are numeric (for example, no NA values)
    if (!is.numeric(xvec)) 
        stop("x must be numeric")
    if (!is.numeric(yvec))
        stop("y must be numeric")

#Check that mean has dimension 1
    if (length(mean) != 1)
        stop("mean must have dimension 1")

#"Clean" data using fcn myWdataclean2="weighted dataclean 2"
    temp1x <- myWdataclean2(xvec, dx)
    temp1y <- myWdataclean2(yvec, dy)

#Redefine x,y as ascending, distinct data values
    x <- temp1x$value
    y <- temp1y$value
#Redefine status dx,dy corresponding to redefined x,y
    dx <- temp1x$dd
    dy <- temp1y$dd
#Define weights wx,wy as the number of data at each distinct value
    wx <- temp1x$weight
    wy <- temp1y$weight
#Set highest of all status 1,0 to 1
    xindex10 <- which(dx != 2)
    yindex10 <- which(dy != 2)
    dx[xindex10[length(xindex10)]] <- 1
    dy[yindex10[length(yindex10)]] <- 1
#Set lowest of all status 1,2 to 1
    xindex12 <- which(dx != 0)
    yindex12 <- which(dy != 0)
    dx[xindex12[1]] <- 1
    dy[yindex12[1]] <- 1
#Put censored,uncensored data into respective vectors 
    xd0  <-  x[dx == 0]
    wxd0 <- wx[dx == 0]
    xd1  <-  x[dx == 1]
    wxd1 <- wx[dx == 1]
    xd2  <-  x[dx == 2]
    wxd2 <- wx[dx == 2]
    yd0  <-  y[dy == 0]
    wyd0 <- wy[dy == 0]
    yd1  <-  y[dy == 1]
    wyd1 <- wy[dy == 1]
    yd2  <-  y[dy == 2]
    wyd2 <- wy[dy == 2]
#Check that there are at least 2 uncensored data for x,y 
    if (length(xd1) <= 1) 
        stop("need more distinct uncensored x obs.")
    if (length(yd1) <= 1) 
        stop("need more distinct uncensored y obs.")
#Store vector lengths 
    nx0 <- length(xd0)
    ny0 <- length(yd0)
    nx1 <- length(xd1)
    ny1 <- length(yd1)
    nx2 <- length(xd2)
    ny2 <- length(yd2)
    nx  <- length(x)
    ny  <- length(y)  

#LR denominator (unconstrained):
   #Doubly-censored case
    if ( (nx0>0) & (nx2>0 ) ) {
      temp2x   <- myWCY(x = x, d = dx, wt = wx)
      logelx00 <- temp2x$logEL
      jumpxu   <- temp2x$jump
   #Adjust xd1, nx1 etc. to reflect repeat data
      xd1<-temp2x$xd1
      nx1<-length(xd1)
      wxd1<-temp2x$wd1
      xd0<-temp2x$xd0
      nx0<-length(xd0)
      wxd0<-temp2x$wd0
      }
    if ( (ny0>0) & (ny2>0 ) ) {
      temp2y   <- myWCY(x = y, d = dy, wt = wy)
      logely00 <- temp2y$logEL
      jumpyu   <- temp2y$jump
   #Adjust yd1, ny1 etc. to reflect repeat data
      yd1<-temp2y$xd1
      ny1<-length(xd1)
      wyd1<-temp2y$wd1
      yd0<-temp2y$xd0
      ny0<-length(xd0)
      wyd0<-temp2y$wd0
      }
   #Right-censored case
    if ( (nx0>0) & (nx2==0 ) ) {
      temp2x   <- myWKM(x = x, d = dx, w = wx)
      logelx00 <- temp2x$logel
      jumpxu   <- temp2x$jump[temp2x$jump>0]
   #Adjust xd1, nx1 etc. to reflect repeat data
      xd1<-temp2x$times[temp2x$jump>0]
      nx1<-length(xd1)
      wxd1<-temp2x$weight[temp2x$jump>0]
      xd0<-temp2x$times[temp2x$jump==0]
      nx0<-length(xd0)
      wxd0<-temp2x$weight[temp2x$jump==0]
      }
    if ( (ny0>0) & (ny2==0 ) ) {
      temp2y   <- myWKM(x = y, d = dy, w = wy)
      logely00 <- temp2y$logel
      jumpyu   <- temp2y$jump[temp2y$jump>0]
   #Adjust yd1, ny1 etc. to reflect repeat data
      yd1<-temp2y$times[temp2y$jump>0]
      ny1<-length(yd1)
      wyd1<-temp2y$weight[temp2y$jump>0]
      yd0<-temp2y$times[temp2y$jump==0]
      ny0<-length(yd0)
      wyd0<-temp2y$weight[temp2y$jump==0]
      } 
   #Left-censored case
    if ( (nx0==0) & (nx2>0 ) ) { 
      dlx <- dx 
      dlx [dlx == 2] <- 0     
      temp2x   <- myWKM(x = x, d = rev(dlx), w = rev(wx))
      logelx00 <- temp2x$logel
      jumpxu   <- rev(temp2x$jump[temp2x$jump>0])
   #Adjust xd1, nx1 etc. to reflect repeat data
      xd1<-temp2x$times[rev(temp2x$jump)>0]
      nx1<-length(xd1)
      wxd1<-rev(temp2x$weight)[rev(temp2x$jump)>0]
      xd2<-temp2x$times[rev(temp2x$jump)==0]
      nx2<-length(xd2)
      wxd2<-rev(temp2x$weight)[rev(temp2x$jump)==0]
      }
    if ( (ny0==0) & (ny2>0 ) ) {
      dly <- dy
      dly[dly == 2] <- 0
      temp2y   <- myWKM(x = y, d = rev(dly), w = rev(wy))
      logely00 <- temp2y$logel
      jumpyu   <- rev(temp2y$jump[temp2y$jump>0])
   #Adjust yd1, ny1 etc. to reflect repeat data
      yd1<-temp2y$times[rev(temp2y$jump)>0]
      ny1<-length(yd1)
      wyd1<-rev(temp2y$weight)[rev(temp2y$jump)>0]
      yd2<-temp2y$times[rev(temp2y$jump)==0]
      ny2<-length(yd2)
      wyd2<-rev(temp2y$weight)[rev(temp2y$jump)==0]
      }
   #Uncensored case
    if ( (nx0==0) & (nx2==0 ) ) { 
      logelx00 <- sum(wxd1 * log(wxd1/sum(wxd1)))
      jumpxu   <- wxd1/sum(wxd1)
         }
    if ( (ny0==0) & (ny2==0 ) ) {
      logely00 <- sum(wyd1 * log(wyd1/sum(wyd1)))
      jumpyu   <- wyd1/sum(wyd1)
         }
   #Calculate likelihood
      logel00 <- logelx00 + logely00
   #Calculate NPMLE
      indic <- matrix(NA,nrow=nx1,ncol=ny1)
      for (i in 1:nx1) {
       for (j in 1:ny1) {
         indic[i,j] <- fun(xd1[i],yd1[j]) } }
      indicmat <- indic - mean
      funNPMLE=as.vector(jumpxu %*% indic %*% jumpyu)

#LR numerator (constrained): 
   #Initialize muvec,nuvec
    muvec <- jumpxu
    nuvec <- jumpyu
   #Initialize kx,ky and kkx,kky
    #For each rt-cens x,y value, kx,ky holds xd1,yd1 index of first
    #uncensored value to its right
    #For each left-censored value, kkx,kky holds xd1,yd1 index of first
    #uncensored value to its left 
     if (nx0>0) {
       kx <- rep(NA, nx0)
        for (i in 1:nx0) {kx[i] <- 1 + nx1 - sum(xd1 > xd0[i])}
      }
     if (nx2>0) {
       kkx <- rep(NA, nx2)
       for (i in 1:nx2) {kkx[i] <- sum(xd1 < xd2[i])}
      }
     if (ny0>0) {
       ky <- rep(NA, ny0)
        for (j in 1:ny0) {ky[j] <- 1 + ny1 - sum(yd1 > yd0[j])}
      }
     if (ny2>0) {
       kky <- rep(NA, ny2)
       for (j in 1:ny2) {kky[j] <- sum(yd1 < yd2[j])}
      }
  #Initialize iteration counter num and log-likelihood-holder logvec
     num <- 1
     logvec <- rep(0,maxit)
  #Repeat EM repeatedly for maxit iterations
       while (num <= maxit) {
  #Initialize weights wxd1new,wyd1new
        wxd1new <- wxd1
        wyd1new <- wyd1
  #Perform Expectation step on uncensored x,y weights
        if (nx0>0) {
          surx <- rev(cumsum(rev(muvec)))
          for (i in 1:nx0) {wxd1new[kx[i]:nx1] <- wxd1new[kx[i]:nx1] +
            wxd0[i] * muvec[kx[i]:nx1]/surx[kx[i]] }
          }
        if (nx2>0) {
          cdfx <- cumsum(muvec)
          for (i in 1:nx2) {wxd1new[1:kkx[i]] <- wxd1new[1:kkx[i]] +
            wxd2[i] * muvec[1:kkx[i]]/cdfx[kkx[i]] }
          }       	
        if (ny0>0) {
          sury <- rev(cumsum(rev(nuvec)))
          for (j in 1:ny0) {wyd1new[ky[j]:ny1] <- wyd1new[ky[j]:ny1] +
            wyd0[j] * nuvec[ky[j]:ny1]/sury[ky[j]] }
           }
        if (ny2>0) {
          cdfy <- cumsum(nuvec)
          for (j in 1:ny2) {wyd1new[1:kky[j]] <- wyd1new[1:kky[j]] +
            wyd2[j] * nuvec[1:kky[j]]/cdfy[kky[j]] }
           }           
  #Perform Maximization step on uncensored x,y jumps
        temp3 <- el2.test.wts(xd1, yd1, wxd1new, wyd1new, muvec,
          nuvec, indicmat, mean)
        muvec <- temp3$jumpu
        nuvec <- temp3$jumpv        
   #Calculate loglikelihood so its convergence can be tracked
        logelx <- sum(wxd1 * log(muvec))
         if (nx0>0) {
          surx <- rev(cumsum(rev(muvec)))
          logelx <- logelx + sum(wxd0 * log(surx[kx]))
         }
         if (nx2>0) {
          cdfx <- cumsum(muvec)
          logelx <- logelx + sum(wxd2 * log(cdfx[kkx]))
         }
       logely <- sum(wyd1 * log(nuvec))
         if (ny0>0) {
          sury <- rev(cumsum(rev(nuvec)))
          logely <- logely + sum(wyd0 * log(sury[ky]))
         }
         if (ny2>0) {
          cdfy <- cumsum(nuvec)
          logely <- logely + sum(wyd2 * log(cdfy[kky]))
         }
       logel <- logelx + logely
       logvec[num]<-logel
       num <- num + 1
       }
    #Store -2log(likelihood ratio) in tval
       tval <- 2 * (logel00 - logel)
    #Calculate constraint
       constraint <- as.vector(muvec %*% indicmat %*% t(nuvec))
   #Return results of the el2.cen.EM function 
       list(xd1=xd1,yd1=yd1,temp3=temp3, mean=mean, funNPMLE=
        funNPMLE, logel00=logel00, logel=logel, "-2LLR"=tval,
        Pval=1 - pchisq(tval, df = 1), logvec=logvec,
        sum_muvec=sum(muvec), sum_nuvec=sum(nuvec),
        constraint=constraint)
}