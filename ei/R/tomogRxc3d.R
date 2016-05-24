# 3-d Tomography Plot for 2x3 Case

tomogRxC3d <- function(formula,data,total=NULL,lci=TRUE,estimates=FALSE,ci=FALSE,level=.95,seed=1234,color=hcl(h=30,c=100,l=60),transparency=.75,light=FALSE,rotate=TRUE){
  ##Run Through RxC Code Once (from Molly's original tomogRxC function)
  #require(grDevices)

  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package rgl is needed for the tomogRxC3d function to work. Please install it.",
      call. = FALSE)
  }
  
  noinfocount <- 0
  form <- formula
  dvname <- terms.formula(formula)[[2]]
  covariate <- NA
  #Make the bounds
  rows <- c(all.names(form)[6:(length(all.names(form)))])
  names <- rows
  cols <- c(all.names(form)[3])
  options(warn=-1)
  bnds <- bounds(form, data=data,rows=rows, column =cols,threshold=0)
  options(warn=0)
  #Totals
  dv <- data[, all.names(form)[3]]
  #Assign other category
  bndsoth <- bnds$bounds[[3]]
  oth <- data[,all.names(form)[length(all.names(form))]]
  
  #Assign x-axis category
  bndsx <- bnds$bounds[[1]]
  xcat <- data[,all.names(form)[6]]
  
  #Assign y-axis category
  bndsy <- bnds$bounds[[2]]
  ycat <- data[,all.names(form)[7]]
  
  #Minimums & Maximums
  minx <- bndsx[,1]
  miny <- bndsy[,1]
  minoth <- bndsoth[,1]
  maxx <- bndsx[,2]
  maxy <- bndsy[,2]
  maxoth <- bndsoth[,2]
  
  #####
  #Starting point when x is at minimum
  ##
  #Holding x at its minimum, what are the bounds on y?
  
  #When x is at its minimum, the new dv and total are:
  newdv <- dv - (minx*xcat)
  newtot <- oth + ycat 
  t <- newdv/newtot
  y <- ycat/newtot
  
  #The new bounds on the y category are:
  
  lby <- cbind(miny, (t - maxoth*oth/newtot)/(y))
  lby[,2] <- ifelse(y==0, 0, lby[,2])
  lowy <- apply(lby,1,max)
  hby <- cbind((t-minoth*oth/newtot)/y,maxy)
  highy <- apply(hby,1,min)
  
  #####
  #Starting point when x is at maximum
  ##
  #Holding x at its maximum, what are the bounds on y?
  
  #The new bounds on x are:
  newtot <- oth + xcat
  newdv <- dv - (miny*ycat)
  x <- xcat/newtot
  t <- newdv/newtot
  lbx <- cbind(minx, (t-maxoth*oth/newtot)/x)
  lbx[,2] <- ifelse(x==0, 0, lbx[,2])
  lowx <- apply(lbx,1,max)
  hbx <- cbind((t-minoth*oth/newtot)/x,maxx)
  highx <- apply(hbx,1,min)
  
  #Graph starting points
  #High starting points
  hstr <- cbind(minx, highy)
  #High ending points
  hend <- cbind(highx, miny)
  
  #Low starting points
  lstr <- cbind(minx, lowy)
  lend <- cbind(lowx, miny)
  
  xl <- paste("Percent", names[1], dvname[2])
  yl <- paste("Percent", names[2], dvname[2])
  zl <- paste("Percent", names[3], dvname[2])
  mn <- paste("3-d Tomography Plot in a 2x3 Table")
  
  #plot(c(0,0), xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i",xlab=xl, ylab=yl, col="white", main=mn)
  
  
  ok <- !is.na(hstr[,2]) & !is.na(hend[,1])
  exp1 <- hstr[ok,2]>=maxy[ok]
  exp2 <- hend[ok,1]>=maxx[ok]
  exp3 <- lstr[ok,2]<=miny[ok]
  exp4 <- lend[ok,1]<=minx[ok]
  hstr <- hstr[ok,]
  hend <- hend[ok,]
  lstr <- lstr[ok,]
  lend <- lend[ok,]
  dv <- dv[ok]
  ycat <- ycat[ok]
  oth <- oth[ok]
  minoth <- minoth[ok]
  xcat <- xcat[ok]
  maxy <- maxy[ok]
  maxx <- maxx[ok]
  
  ##Create Matrices for coordinates
  coords.xaxs<-matrix(data=NA,nrow=dim(hstr)[1],ncol=10)
  coords.yaxs<-matrix(data=NA,nrow=dim(hstr)[1],ncol=10)
  
  for(i in 1:dim(hstr)[1]){
    if((exp1[i] + exp2[i] + exp3[i] + exp4[i])==0){
      xaxs <- c(hstr[i,1],  lstr[i,1],lend[i,1],hend[i,1])
      yaxs <- c(hstr[i,2], lstr[i,2],lend[i,2], hend[i,2])
      side1 <- c(xaxs,xaxs[1])
      side2 <- c(yaxs,yaxs[1])
      c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
      area <- abs(c1-c2)/2
      #polygon(xaxs,yaxs, col=rgb(0,0,0,alpha=0), border="black", lty=1, lwd=.25)
    }
    if((exp1[i]==1) & (exp2[i])==0){
      cut <- (dv[i]-(oth[i])*minoth[i])/xcat[i] - maxy[i]*ycat[i]/xcat[i]
      kink1x <- c(cut)
      kink1y <- c(maxy[i])
    }
    if((exp2[i]==1) & (exp1[i])==0){
      cut <- (dv[i]-(oth[i])*minoth[i])/ycat[i] - maxx[i]*xcat[i]/ycat[i]
      kink1x <- c(maxx[i])
      kink1y <- c(cut)
    }
    if((exp2[i]==1 & exp1[i]==1)){
      cut <- (dv[i]-(oth[i])*minoth[i])/ycat[i] - maxx[i]*xcat[i]/ycat[i]
      cut2 <- (dv[i]-(oth[i])*minoth[i])/xcat[i] - maxy[i]*ycat[i]/xcat[i]
      kink1x <- c(maxx[i], cut2)
      kink1y <- c(cut, maxy[i])
    }
    if((exp3[i]==1) & (exp4[i])==0){
      cut <- (dv[i]-(oth[i])*maxoth[i])/xcat[i] - miny[i]*ycat[i]/xcat[i]
      kink2x <- c(cut)
      kink2y <- c(miny[i])
    }
    if((exp4[i]==1) & (exp3[i])==0){
      cut <- (dv[i]-(oth[i])*maxoth[i])/ycat[i] - minx[i]*xcat[i]/ycat[i]
      kink2x <- c(minx[i])
      kink2y <- c(cut)
    }
    if((exp3[i]==1 & exp4[i]==1)){
      cut <- (dv[i]-(oth[i])*maxoth[i])/ycat[i] - minx[i]*xcat[i]/ycat[i]
      cut2 <- (dv[i]-(oth[i])*maxoth[i])/xcat[i] - miny[i]*ycat[i]/xcat[i]
      kink2x <- c(minx[i], cut2)
      kink2y <- c(cut, miny[i])
    }
    if((exp3[i] + exp4[i])==0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i])!=0){
      xaxs <- c(hstr[i,1],  lstr[i,1],lend[i,1],hend[i,1], kink1x)
      xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
      yaxs <- c(hstr[i,2], lstr[i,2],lend[i,2], hend[i,2], kink1y)
      yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
      side1 <- c(xaxs,xaxs[1])
      side2 <- c(yaxs,yaxs[1])
      c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
      area <- abs(c1-c2)/2
      #polygon(xaxs,yaxs, col=rgb(0,0,0,alpha=0), border="black", lty=1, lwd=.25)
    }
    if((exp1[i] + exp2[i])==0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i])!=0){
      xaxs <- c(hstr[i,1],  lstr[i,1],kink2x,lend[i,1],hend[i,1])
      xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
      yaxs <- c(hstr[i,2], lstr[i,2],kink2y,lend[i,2], hend[i,2])
      yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
      side1 <- c(xaxs,xaxs[1])
      side2 <- c(yaxs,yaxs[1])
      c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
      area <- abs(c1-c2)/2
      #polygon(xaxs,yaxs, col=rgb(0,0,0,alpha=0), border="black", lty=1, lwd=.25)
    }
    if((exp1[i] + exp2[i])!=0 & (exp3[i] + exp4[i])!=0){
      xaxs <- c(hstr[i,1],  lstr[i,1],kink2x,lend[i,1],hend[i,1], kink1x)
      xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
      yaxs <- c(hstr[i,2], lstr[i,2],kink2y,lend[i,2], hend[i,2], kink1y)
      yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
      side1 <- c(xaxs,xaxs[1])
      side2 <- c(yaxs,yaxs[1])
      c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
      area <- abs(c1-c2)/2
      #polygon(xaxs,yaxs, col=rgb(0,0,0,alpha=0), border="black", lty=1, lwd=.225)
    }
    for(j in 1:length(xaxs)){if(xaxs[j]<1e-5){xaxs[j]<-0}}
    for(k in 1:length(yaxs)){if(yaxs[k]<1e-5){yaxs[k]<-0}}
    for(l in 1:length(xaxs)){if(xaxs[l]>.99999){xaxs[l]<-1}}
    for(m in 1:length(yaxs)){if(yaxs[m]>.99999){yaxs[m]<-1}}
    axs<-rbind(xaxs,yaxs)
    coords.axs<-unique(axs,MARGIN=2,incomparables=F)
    coords.axs<-if(length(coords.axs[1,])!=10){cbind(coords.axs,matrix(data=NA,nrow=2,ncol=(10-length(coords.axs[1,]))))}else{coords.axs}
    coords.xaxs[i,]<-coords.axs[1,]
    coords.yaxs[i,]<-coords.axs[2,]
  }
  
  ##Record Coordinate
  w.coord<-coords.xaxs[,colSums(is.na(coords.xaxs))<nrow(coords.xaxs)]
  b.coord<-coords.yaxs[,colSums(is.na(coords.yaxs))<nrow(coords.yaxs)]
  
  
  #######################################################################
  
  ##THIS WILL BE A FN WHERE YOU SPECIFY NAME OF TURNOUT, TOTAL, CATEGORIES
  
  ##Combine Coordinates for 3d Graph
  w<-matrix(NA,nrow=dim(hstr)[1],ncol=dim(w.coord)[2])
  b<-matrix(NA,nrow=dim(hstr)[1],ncol=dim(w.coord)[2])
  h<-matrix(NA,nrow=dim(hstr)[1],ncol=dim(w.coord)[2])
  
  #Extract Data needed for accounting Identity
  N <- data[,all.names(form)[3]] + data[,all.names(form)[4]]
  t <- data[,all.names(form)[3]]/N
  xb <- data[,all.names(form)[7]]/N
  xw <- data[,all.names(form)[6]] /N
  
  #Using the Accounting Identity from King (1997), calculate coordinates for remaining category
  ##T = B_b*X_b + B_w*X_w + B_h*X_h = B_b*X_b + B_w*X_w + B_h*(1-X_b-X_h)
  coords<-matrix(NA, nrow=dim(w.coord)[2], ncol=3)
  for(i in 1: dim(hstr)[1]){
    coords<-cbind(w.coord[i,],b.coord[i,],rep(NA, nrow(coords)))  
    for(j in 1:nrow(coords)){
      #Account for cases where xh=0
      if(1-xb[i]-xw[i] < 1e-5){coords[j,3]<-0}else{
        coords[j,3]<-(t[i] - coords[j,2]*xb[i] - coords[j,1]*xw[i])/(1-xb[i]-xw[i])}
      ##Make Sure that it records small numbers as 0
      if(coords[j,3]<1e-5 & is.na(coords[j,3])==F){coords[j,3]<-0}
      if(coords[j,3]>.99999 & is.na(coords[j,3])==F){coords[j,3]<-1}
      ##Account for cases where xb=0 and xw
      if(xb[i]==0){coords[j,2]<-0}
      if(xw[i]==0){coords[j,1]<-0}
    }
    w[i,]<-coords[,1]
    b[i,]<-coords[,2]
    h[i,]<-coords[,3]
  }
  
  ##Run RxC Model
  set.seed(seed)
  dbuf <- ei(formula=formula, data=data)
  
  ##Calculate Center of Ellipsoid: mean of estimates
  means<-apply(dbuf$draws$Beta, MARGIN=2, FUN=mean) #Mean of all simulations for each precinct and beta
  means.w <- means[seq(1, length(means), 6)]
  means.b <- means[seq(2, length(means), 6)]
  means.h <- means[seq(3, length(means), 6)]
  
  ##Calculate Variance of Ellipsoid: average variance of simulations
  ses.w <- sd(means.w)
  ses.b <- sd(means.b)
  ses.h <- sd(means.h)
  
  ##Calculate Covariance of Ellipsoid (Rho): average covariance of simulations
  covs.bw<-cov(means.w,means.b)
  covs.wh<-cov(means.w,means.h)
  covs.bh<-cov(means.b,means.h)
  #Create Rho Matrix
  #Rho if covariance divided by multiplied variance
  rho<-matrix(c(1,covs.bw/(ses.b*ses.w),covs.bh/(ses.b*ses.h),covs.bw/(ses.b*ses.w),1,covs.wh/(ses.h*ses.w)
                ,covs.bh/(ses.b*ses.h),covs.wh/(ses.h*ses.w),1),nrow=3)
  
  
  ######################################
  
  ##CALCULATE AREA OF PLANES for LoCoI
  area <- rep(NA, nrow(w))
  for(i in 1:nrow(w)){
    #2-sided Line
    if(length(unique(w[i,]))==1 | length(unique(b[i,]))==1 | length(unique(h[i,]))==1){
      if(length(na.omit(w[i,]))>2){area[i]<-(sqrt(abs(h[i,1]-h[i,2])^2 + abs(b[i,1]-b[i,2])^2))*.01} #Multiply by epsilon to weight it appropriately small
      if(length(na.omit(b[i,]))>2){area[i]<-(sqrt(abs(h[i,1]-h[i,2])^2 + abs(w[i,1]-w[i,2])^2))*.01}
      if(length(na.omit(h[i,]))>2){area[i]<-(sqrt(abs(w[i,1]-w[i,2])^2 + abs(b[i,1]-b[i,2])^2))*.01}
    }
    # All other Planes
    if(length(unique(w[i,]))!=1 & length(unique(b[i,]))!=1 & length(unique(h[i,]))!=1){
      # Get 3d coordinates
      p <- cbind(b[i,],w[i,],h[i,])
      p <- na.omit(p)
      
      # Get equation of plane
      fit <- lm(h[i,] ~ b[i,] + w[i,])
      coefs <- coef(fit)
      a.eq <- coefs[2]
      b.eq <- coefs[3]
      c.eq <- -1
      d.eq <- coefs["(Intercept)"]
      
      # Get normal vector
      n <- matrix(c(a.eq,b.eq,c.eq), nrow=1, ncol=3)
      
      # Calculate Basis Vector on Plane
      u1 <- matrix(c( (p[2,1]-p[1,1]), (p[2,2]-p[1,2]), (p[2,3]-p[1,3])), nrow=1, ncol=3)
      
      # Get Second Basis Vector as cross-product of n and u1
      CrossProduct3D <- function(x, y, i=1:3) {
        # Project inputs into 3D, since the cross product only makes sense in 3D.
        To3D <- function(x) head(c(x, rep(0, 3)), 3)
        x <- To3D(x)
        y <- To3D(y)
        # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
        # so on).  Index3D() lets us do that using R's convention of 1-based (rather
        # than 0-based) arrays.
        Index3D <- function(i) (i - 1) %% 3 + 1
        # The i'th component of the cross product is:
        # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
        # as long as we treat the indices cyclically.
        return (x[Index3D(i + 1)] * y[Index3D(i + 2)] -
                  x[Index3D(i + 2)] * y[Index3D(i + 1)])
      }
      
      u2 <- matrix(CrossProduct3D(n,u1), nrow=1, ncol=3)
      
      # Normalize the basis vectors
      e1 <- u1/sqrt(u1[1,1]^2 + u1[1,2]^2 + u1[1,3]^2)
      e2 <- u2/sqrt(u2[1,1]^2 + u2[1,2]^2 + u2[1,3]^2)
      
      # Get the projection matrix from basis vectors
      pm <- rbind(e1, e2)
      
      # Get the 2d coordinates from projection matrix
      p.2d <- t(pm%*%t(p))
      
      # Create matrix of points in counter-clockwise order (need to systematize)
      ccorder <- rev(chull(p.2d)) #Convex hull gives clockwise order then reverse
      parea <- p.2d[order(ccorder),]
      
      # Add first row to bottom of matrix
      parea <- rbind(parea, parea[1,])
      
      # Calcluate area as 1/2 of determinant (need to systematize)
      plus.area <- 0
      for(j in 2:dim(parea)[1]){
        x <- parea[j-1,1]*parea[j,2]
        plus.area <- plus.area + x
      }
      minus.area <- 0
      for(k in 2:dim(parea)[1]){
        x <- parea[k-1,2]*parea[k,1]
        minus.area <- minus.area + x
      }
      area[i] <- 1/2*(plus.area-minus.area)
    }
  }
  
  scale<-(area-min(area))/(max(area)-min(area))*100
  
  #################################
  ## Plot Based on Plot Type
  
  ##Tomography Plot - with/without point estimates
  if(ci==F){
    rgl::plot3d(b,w,h, type="s", col="black", size=0.3, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), xlab=yl, ylab=xl, zlab=zl)
    ##Add the plane for each set of points
    for(i in 1:dim(hstr)[1]){  
      #Two point plane (line)
      if(length(unique(w[i,]))==1 | length(unique(b[i,]))==1 | length(unique(h[i,]))==1){
        if(lci==T){rgl::lines3d(x = b[i,], y = w[i,], z=h[i,], col=hcl(h=30,c=100,l=scale[i],alpha=1))}
        else{rgl::lines3d(x = b[i,], y = w[i,], z=h[i,], col=color)}
      }else{
        #All Others
        fit <- lm(h[i,] ~ b[i,] + w[i,])
        coefs <- coef(fit)
        a <- coefs[2]
        c <- coefs[3]
        d <- -1
        e <- coefs["(Intercept)"]
        if(lci==T){rgl::planes3d(a, c, d, e, alpha=transparency, col=hcl(h=30,c=100,l=scale[i]), lit=light)}
        else{rgl::planes3d(a, c, d, e, alpha=transparency, col=color, lit=light)}
      }
    }
    # With point estimates - mean of average simulations across precincts
    if(estimates==T){rgl::spheres3d(means.b,means.w, means.h, col="black", radius=.03, add=T)}
  }
  
  ##With Confidence Elipse
  if(ci==T){  
    rgl::plot3d(rgl::ellipse3d(x=rho , scale= c(ses.b,ses.w,ses.h), centre= c(mean(means.b),mean(means.w), mean(means.h))
                     , alpha=.8, level=level) ,color=hcl(h=10,c=60,l=20), xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), xlab=yl, ylab=xl, zlab=zl, lit=light)
    rgl::points3d(mean(means.b),mean(means.w), mean(means.h), col="black", size=10, add=T) #Mean of average simulations across precincts
    rgl::points3d(b,w,h, col="black", size=1, add=T)
    ##Add the plane for each set of points
    for(i in 1:dim(hstr)[1]){  
      #Two point plane (line)
      if(length(unique(w[i,]))==1 | length(unique(b[i,]))==1 | length(unique(h[i,]))==1){
        if(lci==T){rgl::lines3d(x = b[i,], y = w[i,], z=h[i,], col=hcl(h=30,c=100,l=scale[i],alpha=1))}
        else{rgl::lines3d(x = b[i,], y = w[i,], z=h[i,], col=hcl(h=30,c=100,l=60,alpha=1))}
      }else{
        #All Others
        fit <- lm(h[i,] ~ b[i,] + w[i,])
        coefs <- coef(fit)
        a <- coefs[2]
        c <- coefs[3]
        d <- -1
        e <- coefs["(Intercept)"]
        if(lci==T){rgl::planes3d(a, c, d, e, alpha=.25, col=hcl(h=30,c=100,l=scale[i]), lit=light)}
        else{rgl::planes3d(a, c, d, e, alpha=.25, col=hcl(h=30,c=100,l=60), lit=light)}
      }
    }
    # With point estimates - mean of average simulations across precincts
    if(estimates==T){rgl::spheres3d(means.b,means.w, means.h, col="black", radius=.03, add=T)}
  }
  
  if(rotate==T){rgl::play3d( rgl::spin3d(rpm=2.5), duration=20)}
}