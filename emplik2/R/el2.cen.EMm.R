el2.cen.EMm<-function(x, dx, y, dy, p, H, xc=1:length(x),
  yc=1:length(y), mean, maxit=10) {
#x,y pairs can be any combination of uncensored, left-cens, right-cens
#Store data as vectors 
    xvec <- as.vector(x)
    yvec <- as.vector(y)
    dx   <- as.vector(dx)
    dy   <- as.vector(dy)
    mean <- as.vector(mean)
    xc   <- as.vector(xc)
    yc   <- as.vector(yc)
    nx   <- length(xvec)
    ny   <- length(yvec)
#Check that status and observations have same length in x,y
    if (length(dx) != nx) 
        stop("length of x and dx must agree")
    if (length(dy) != ny) 
        stop("length of y and dy must agree")
    if (length(xc) != nx) 
        stop("length of xc and dx must agree")
    if (length(dy) != ny) 
        stop("length of yc and dy must agree")
#Check that there are enough data in x,y
    if (nx <= length(mean)) 
        stop("need more observations than length of mean in x")
    if (ny <= 2*length(mean) + 1) 
        stop("need more observations than length of mean in y")
#Check that status are only 0,1,2
    if (any((dx != 0) & (dx != 1) & (dx != 2)))  
        stop("dx must be 0(right-censored) or 1(uncensored) or 2(left-
         censored")
    if (any((dy != 0) & (dy != 1) & (dy != 2)))  
        stop("dy must be 0(right-censored) or 1(uncensored) or 2(left-
         censored")    
#Check that xvec,yvec are numeric
    if (!is.numeric(xvec)) 
        stop("x must be numeric")
    if (!is.numeric(yvec))
        stop("y must be numeric")
    if (!is.numeric(mean))
        stop ("mean must be numeric")
#Check that H dimensions are consistent with lengths of x,y
    if (dim(H)[1] != length(x) )
        stop("dim(H)[1] must equal length(x)")
    if (dim(H)[2] != p*length(y) )
        stop("dim(H)[2] must equal p*length(y)")
#"Clean" data using code from Wdataclean5 ("Weighted Dataclean5")
    #Clean x data
    niceorderx <- order(x, -dx)
    x   <-  x[niceorderx]
    dx  <- dx[niceorderx]
    wx  <- wx[niceorderx]
    xc  <- xc[niceorderx]
    t1  <- x[-1]  != x[-nx]
    t2  <- dx[-1] != dx[-nx]
    t3  <- xc[-1] != xc[-nx]
    t <- t1 | t2 | t3
    ind <- c(which(t | is.na(t)), nx)
    csumwx <- cumsum(wx)
    x  <- x[ind]
    dx <- dx[ind]
    wx <- diff(c(0, csumwx[ind]))
    H   <- as.matrix(H[niceorderx, ][ind,])
    #Clean y data
    niceordery <- order(y, -dy)  
    y   <-  y[niceordery]
    dy  <- dy[niceordery]
    wy  <- wy[niceordery]
    yc  <- yc[niceordery]
    t1  <- y[-1]  != y[-ny]
    t2  <- dy[-1] != dy[-ny]
    t3  <- yc[-1] != yc[-ny]
    t <- t1 | t2 | t3
    ind <- c(which(t | is.na(t)), ny)
    csumwy <- cumsum(wy)
    y  <- y[ind]
    dy <- dy[ind]
    wy <- diff(c(0, csumwy[ind]))        
    for (k in 1:p) {      
        H[, ((k-1)*ny+1):(k*ny)] <- 
          H[, ((k-1)*ny+1):(k*ny)][,niceordery][,ind] }     
#For data with status 0 or 1, set status of highest datum to 1
    xindex10 <- which(dx != 2)
    yindex10 <- which(dy != 2)
    dx[xindex10[length(xindex10)]] <- 1
    dy[yindex10[length(yindex10)]] <- 1
# For data with status 2 or 1, set status of lowest datum to 1
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
      }
    if ( (ny0>0) & (ny2>0 ) ) {
      temp2y   <- myWCY(x = y, d = dy, wt = wy)
      logely00 <- temp2y$logEL
      jumpyu   <- temp2y$jump
      }
   #Right-censored case
    if ( (nx0>0) & (nx2==0 ) ) {
      temp2x   <- myWKM(x = x, d = dx, w = wx)
      logelx00 <- temp2x$logel
      jumpxu   <- temp2x$jump[temp2x$jump>0]
      }
    if ( (ny0>0) & (ny2==0 ) ) {
      temp2y   <- myWKM(x = y, d = dy, w = wy)
      logely00 <- temp2y$logel
      jumpyu   <- temp2y$jump[temp2y$jump>0]
      } 
   #Left-censored case
    if ( (nx0==0) & (nx2>0 ) ) { 
      dlx <- dx 
      dlx [dlx == 2] <- 0     
      temp2x   <- myWKM(x = x, d = rev(dlx), w = rev(wx))
      logelx00 <- temp2x$logel
      jumpxu   <- rev(temp2x$jump[temp2x$jump>0])
      }
    if ( (ny0==0) & (ny2>0 ) ) {
      dly <- dy
      dly[dly == 2] <- 0
      temp2y   <- myWKM(x = y, d = rev(dly), w = rev(wy))
      logely00 <- temp2y$logel
      jumpyu   <- rev(temp2y$jump[temp2y$jump>0])
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
      NPMLE=rep(NA,p)
      H1 <- as.matrix(H[which(dx==1),])
        for (j in 1:p) {
        H2 <- H1[,((j-1)*ny+1):(j*ny)]
        H2 <- H2[,which(dy==1)]
        NPMLE[j]=jumpxu%*%H2%*%jumpyu
        }
      #LR numerator (constrained):
   #Create mean-centered H, denoted as Hmc
   Hmc <- matrix(NA, nrow=nx, ncol=p*ny)   
   for (k in 1:p) {
      M <- matrix(mean[k], nrow=nx, ncol=ny)
      Hmc[,((k-1)*ny+1):(k*ny)] <-
         H[,((k-1)*ny+1):(k*ny)] - M }

   #Create uncensored version of Hmc, denoted as Hu 
    whichdx <- which(dx == 1)
    whichdy <- which(dy == 1)
    Hu <- matrix(NA, nrow=nx1,ncol=p*ny1)
      for (k in 1:p) {
        Hu[,((k-1)*ny1+1):(k*ny1)] <-
         Hmc[,((k-1)*ny+1):(k*ny)][whichdx,whichdy] }    

   #Calculate Hmu and Hnu matrices (mean-centered)
    Hmu <- matrix(NA,nrow=p, ncol=ny1*nx1)
    Hnu <- matrix(NA,nrow=p, ncol=ny1*nx1) 
    for (i in 1:p) {
      for (k in 1:nx1) {
        Hmu[i, ((k-1)*ny1+1):(k*ny1)] <- 
          Hu[k,((i-1)*ny1+1):(i*ny1)] } }
    for (i in 1:p) {
      for (k in 1:ny1) {
        Hnu[i,((k-1)*nx1+1):(k*nx1)] <- Hu[(1:nx1),(i-1)*ny1+k]} }
    
   #Initialize muvec,nuvec using unconstrained jumps
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
        temp3 <- el2.test.wtm(xd1,yd1,wxd1new, wyd1new, muvec,
          nuvec, Hu, Hmu, Hnu, p, mean) 
        muvec <- temp3$muvec1
        nuvec <- temp3$nuvec1        
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
       logvec[num]<- logel
       num <- num + 1
       }
    #Store -2log(likelihood ratio) in tval
       tval <- 2 * (logel00 - logel)
    #Return results of the el2.cen.EM function 
       list(xd1=xd1,yd1=yd1,temp3=temp3, NPMLE=
        NPMLE, logel00=logel00, logel=logel, "-2LLR"=tval,
        Pval=1 - pchisq(tval, df = length(mean)), logvec=logvec,
        sum_muvec=sum(muvec), sum_nuvec=sum(nuvec) )
}
