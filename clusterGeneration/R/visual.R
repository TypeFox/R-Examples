# v1.2.5
#  (1) fixed a bug in the function 'plot2DProjection'. The 'lim' should be
#    'ylim' 
#
# plot separation plane a'(y-eta)=0, i.e. a'y-b=0, where b=a'eta.
# y1 -- cluster 1. Rows correspond to observations. 
#                  Columns correspond to variables
# y2 -- cluster 2. Rows correspond to observations. 
#                  Columns correspond to variables
# projDir -- projection direction
# bw -- the smoothing bandwidth to be used by the R function density()
# alpha -- tuning parameter
# xlab, ylab -- labels of x and y coordinates
# title -- title of the plot
plot1DProjection<-function(y1, y2, projDir, 
  sepValMethod=c("normal", "quantile"), bw="nrd0", 
  xlim=NULL, ylim=NULL, 
  xlab="1-D projected clusters", ylab="density estimates", 
  title="1-D Projected Clusters and their density estimates",
  font=2, font.lab=2, cex=1.2, cex.lab=1.2, cex.main=1.5,
  lwd=4, lty1=1, lty2=2, pch1=18, pch2=19, col1=2, col2=4, type="l",
  alpha=0.05, eps=1.0e-10, quiet=TRUE)
{ 
  sepValMethod<-match.arg(sepValMethod, choices=c("normal", "quantile"))
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be between (0, 0.01]!\n")
  }

  x1<-y1%*%projDir
  x2<-y2%*%projDir
  x<-c(x1, x2)
  n1<-length(x1)
  n2<-length(x2)
  n<-n1+n2

  m1<-mean(x1, na.rm=TRUE)
  m2<-mean(x2, na.rm=TRUE)
  if(m1>m2)
  { 
    projDir<- -projDir
    x1<- -x1
    x2<- -x2
    if(!quiet)
    { cat("Warning: Projected cluster 1 is on the right-hand side of projected cluster 2!\n")
     cat("'projDir' is replaced by '-projDir'!\n")
    }
  } 

  x<-c(x1, x2)
  n1<-length(x1)
  n2<-length(x2)

  xx1<-x1
  xx2<-x2
  
  if(sepValMethod=="quantile")
  { if(n1>1) 
    { L1<-quantile(xx1, prob=alpha/2, na.rm=TRUE)
      U1<-quantile(xx1, prob=1-alpha/2, na.rm=TRUE)
    } 
    else 
    { L1<-xx1
      U1<-xx1
    }
    if(n2>1)
    { L2<-quantile(xx2, prob=alpha/2, na.rm=TRUE)
      U2<-quantile(xx2, prob=1-alpha/2, na.rm=TRUE)
    } 
    else 
    { L2<-xx2
      U2<-xx2
    }
  } 
  else 
  { za<-qnorm(1-alpha/2)
    if(n1>1)
    { m1<-mean(xx1, na.rm=TRUE)
      sd1<-sd(c(xx1), na.rm=TRUE)
    } 
    else 
    { m1<-xx1
      sd1<-0
    }
    if(n2>1)
    { m2<-mean(xx2, na.rm=TRUE)
      sd2<-sd(c(xx2), na.rm=TRUE)
    } 
    else 
    { m2<-xx2
      sd2<-0
    }
    L1<-m1-za*sd1
    U1<-m1+za*sd1
    L2<-m2-za*sd2
    U2<-m2+za*sd2
  }

  numer<-(L2-U1)
  denom<-(U2-L1)
  if(abs(denom)<eps) # denominator equal to zero
  { sepVal<- -1 } 
  else { sepVal<-numer/denom }
  
  xlab<-paste(xlab," (sepVal=",round(sepVal,2),", alpha=", round(alpha,2), ")",sep="")

  # obtain density estimation
  if(n1>1)
  { tmp1<-density(xx1, bw=bw)
    tmpx1<-tmp1$x
    tmpy1<-tmp1$y
  } 
  else { tmpx1<-xx1; tmpy1<-1; }
  if(n2>1)
  { tmp2<-density(xx2, bw=bw)
    tmpx2<-tmp2$x
    tmpy2<-tmp2$y
  } 
  else  
  { tmpx2<-xx2
    tmpy2<-1; 
  }
  nn1<-length(tmpx1)
  nn2<-length(tmpx2)

  if(is.null(xlim))
  { xlim<-range(c(xx1, xx2, tmpx1, tmpx2)) }
  if(is.null(ylim))
  { ylim<-range(c(tmpy1, tmpy2, tmpy1, tmpy2)) }


  plotCluster(nn1, nn2, tmpx1, tmpy1, tmpx2, tmpy2, 
              xlim, ylim, xlab, ylab, title, 
              font, font.lab, cex, cex.lab, cex.main, lwd, 
              lty1, lty2, pch1, pch2, col1, col2, type)
 
  plotTickLabel(xx1, xx2, L1,U1,L2,U2, 1, font,
                lwd, lty1, lty2, col1, col2) #ticks and labels
  invisible(list(sepVal=sepVal, projDir=projDir))
}

# function called by plot1DProjection and plot2DProjection
plotCluster<-function(nn1, nn2, tmpx1, tmpy1, tmpx2, tmpy2, 
                      xlim, ylim, xlab, ylab, title, 
                      font, font.lab, cex, cex.lab, cex.main, lwd, 
                      lty1, lty2, pch1, pch2, col1, col2, type="l")
{
  if(nn1>1 && nn2>1)
  { plot(tmpx1, tmpy1, type=type, 
         xlim=xlim, ylim=ylim, 
         xlab=xlab, ylab=ylab,
         col=col1, font=font, cex=cex,
         font.lab=font.lab, cex.lab=cex.lab, lwd=lwd, lty=lty1, pch=pch1) 
    title(title, cex.main=cex.main)
    points(tmpx2, tmpy2, type=type, col=col2, lwd=lwd, lty=lty2, pch=pch2)
  } 
  else if (nn1>1 && nn2==1) 
  { plot(tmpx1, tmpy1, type=type, 
         xlim=xlim, ylim=ylim, 
         xlab=xlab, ylab=ylab,
         col=col1, font=font, cex=cex,
         font.lab=font.lab, cex.lab=cex.lab, lwd=lwd, lty=lty1, pch=pch1) 
    title(title, cex.main=cex.main)
    points(tmpx2, tmpy2, type="p", col=col2, lty=lty2, pch=pch2) 
  }  
  else if (nn1==1 && nn2>1) 
  { plot(tmpx2, tmpy2, type=type, 
         xlim=xlim, ylim=ylim, 
         xlab=xlab, ylab=ylab,
         col=col2, font=font, cex=cex,
         font.lab=font.lab, cex.lab=cex.lab, lwd=lwd, lty=lty1, pch=pch1) 
    title(title, cex.main=cex.main)
    points(tmpx1, tmpy1, type="p", col=col1, lty=lty2, pch=pch2) 
  } 
  else 
  { plot(x=c(tmpx1), y=c(tmpy1), type="p", 
         xlim=xlim, ylim=ylim, 
         xlab=xlab, ylab=ylab,
         col=col1, font=font, cex=cex,
         font.lab=font.lab, cex.lab=cex.lab, lty=lty1, pch=pch1) 
    title(title, cex.main=cex.main)
    points(x=c(tmpx2), y=c(tmpy2), type="p", col=col2, lty=lty2, pch=pch2) 
  }
}

# Plot ticks and labels along rotated 1st axis
# x1, x2 -- projected cluster 1 and 2 along the 1st axis
# L1, U1 -- lower and upper 1-alpha/2 percentile of x1
# L2, U2 -- lower and upper 1-alpha/2 percentile of x2
plotTickLabel<-function(x1, x2, L1, U1, L2, U2, axis=1, font=2, lwd=4,
                        lty1=1, lty2=2, col1=2, col2=4)
{ axis(side=axis,at=x1, labels=FALSE, tick=TRUE,col=col1)
  axis(side=axis,at=x2, labels=FALSE, tick=TRUE,col=col2)
  par(mgp=c(3,2,0))
  axis(side=axis,at=c(L1, U1), labels=c("L1", "U1"), tick=TRUE,col=col1,
    lwd=lwd,font=font, lty=lty1)
  axis(side=axis,at=c(L2, U2), labels=c("L2", "U2"), tick=TRUE,col=col2,
    lwd=lwd,font=font, lty=lty2)
  par(mgp=c(3,1,0))
}

# plot 2-D projection
# sepValMethod -- quantile or normal 
#   quantile means use upper and lower alpha quantiles of univariate projections
#   normal means use mean +/- z(alpha)*sd for mean,sd of projections
# y1 -- cluster 1
# y2 -- cluster 2
# alpha -- tuning parameter
# xlab, ylab -- labels of x and y coordinates
# title -- title of the plot
plot2DProjection<-function(y1, y2, projDir, 
  sepValMethod=c("normal", "quantile"), 
  iniProjDirMethod=c("SL", "naive"), 
  projDirMethod=c("newton", "fixedpoint"), 
  xlim=NULL, ylim=NULL, 
  xlab="1st projection direction", 
  ylab="2nd projection direction", 
  title="Scatter plot of 2-D Projected Clusters",
  font=2, font.lab=2, cex=1.2, cex.lab=1, cex.main=1.5,
  lwd=4, lty1=1, lty2=2, pch1=18, pch2=19, col1=2, col2=4, 
  alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  sepValMethod<-match.arg(sepValMethod, choices=c("normal", "quantile"))
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<=0 || !is.integer(ITMAX))
  {
    stop("The maximum iteration number allowed 'ITMAX' should be a positive integer!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  # obtain first two rotated coordinates
  tmp<-getRotateData(y1, y2, projDir, iniProjDirMethod, projDirMethod,
                     alpha, ITMAX, eps, quiet) 
  Q2<-tmp$Q2 # rotation matrix
  cl<-tmp$cl
  x<-tmp$x

  x1<-x[cl==1,1]
  x2<-x[cl==2,1]
  tmpy1<-x[cl==1,2]
  tmpy2<-x[cl==2,2]

  mx1<-mean(x1, na.rm=TRUE)
  mx2<-mean(x2, na.rm=TRUE)
  if(mx1>mx2)
  { 
    x1<- -x1
    x2<- -x2
    Q2[,1]<- -Q2[,1]
    if(!quiet)
    { cat("Warning: Projected cluster 1 is on the right-hand side of projected cluster 2 along the 1st projection direction!\n")
      cat("'projDir' is replaced by '-projDir'!\n")
    }
  } 
  tmpx<-c(x1, x2)

  my1<-mean(tmpy1, na.rm=TRUE)
  my2<-mean(tmpy2, na.rm=TRUE)
  if(my1>my2)
  { 
    tmpy1<- -tmpy1
    tmpy2<- -tmpy2
    Q2[,1]<- -Q2[,1]
    if(!quiet)
    { cat("Warning: Projected cluster 1 is on the right-hand side of projected cluster 2 along the 1st projection direction!\n")
     cat("'projDir' is replaced by '-projDir'!\n")
    }
  } 
  tmpy<-c(tmpy1, tmpy2)

  xx1<-x1
  xx2<-x2
  yy1<-tmpy1
  yy2<-tmpy2
  n1<-length(xx1)
  n2<-length(xx2)

  if(sepValMethod=="quantile")
  { if(n1>1) 
    { Lx1<-quantile(xx1, prob=alpha/2, na.rm=TRUE)
      Ux1<-quantile(xx1, prob=1-alpha/2, na.rm=TRUE)
      Ly1<-quantile(yy1, prob=alpha/2, na.rm=TRUE)
      Uy1<-quantile(yy1, prob=1-alpha/2, na.rm=TRUE)
    } 
    else 
    { Lx1<-xx1
      Ux1<-xx1
      Ly1<-yy1
      Uy1<-yy1
    }
    if(n2>1)
    { Lx2<-quantile(xx2, prob=alpha/2, na.rm=TRUE)
      Ux2<-quantile(xx2, prob=1-alpha/2, na.rm=TRUE)
      Ly2<-quantile(yy2, prob=alpha/2, na.rm=TRUE)
      Uy2<-quantile(yy2, prob=1-alpha/2, na.rm=TRUE)
    } 
    else 
    { Lx2<-xx2
      Ux2<-xx2
      Ly2<-yy2
      Uy2<-yy2
    }
  } 
  else 
  { za<-qnorm(1-alpha/2)
    if(n1>1)
    { mx1<-mean(xx1, na.rm=TRUE)
      sdx1<-sd(c(xx1), na.rm=TRUE)
      my1<-mean(yy1, na.rm=TRUE)
      sdy1<-sd(c(yy1), na.rm=TRUE)
    } 
    else 
    { mx1<-xx1
      sdx1<-0
      my1<-yy1
      sdy1<-0
    }
    if(n2>1)
    { mx2<-mean(xx2, na.rm=TRUE)
      sdx2<-sd(c(xx2), na.rm=TRUE)
      my2<-mean(yy2, na.rm=TRUE)
      sdy2<-sd(c(yy2), na.rm=TRUE)
    } 
    else 
    { mx2<-xx2
      sdx2<-0
      my2<-yy2
      sdy2<-0
    }
    Lx1<-mx1-za*sdx1
    Ux1<-mx1+za*sdx1
    Lx2<-mx2-za*sdx2
    Ux2<-mx2+za*sdx2
    Ly1<-my1-za*sdy1
    Uy1<-my1+za*sdy1
    Ly2<-my2-za*sdy2
    Uy2<-my2+za*sdy2
  }

  numerx<-(Lx2-Ux1)
  denomx<-(Ux2-Lx1)
  numery<-(Ly2-Uy1)
  denomy<-(Uy2-Ly1)
  if(abs(denomx)<eps) # denominator equal to zero
  { sepValx<- -1 } 
  else { 
    sepValx<-numerx/denomx 
  }
  if(abs(denomy)<eps) # denominator equal to zero
  { sepValy<- -1 } 
  else { sepValy<-numery/denomy }
  
  xlab<-paste(xlab," (sepValx=",round(sepValx,2),", alpha=", round(alpha,2), ")",sep="")
  ylab<-paste(ylab," (sepValy=",round(sepValy,2),", alpha=", round(alpha,2), ")",sep="")

  if(missing(xlim)) { xlim<-range(tmpx) } 
  else { xlim<-xlim }
  if(missing(ylim)) { ylim<-range(tmpy) } 
  else { ylim<-ylim }

  xlim[1]<-min(xlim[1], Lx1, Lx2)
  xlim[2]<-max(xlim[2], Ux1, Ux2)
  ylim[1]<-min(ylim[1], Ly1, Ly2)
  ylim[2]<-max(ylim[2], Uy1, Uy2)

  plotCluster(n1, n2, xx1, yy1, xx2, yy2, 
              xlim, ylim, xlab, ylab, title, 
              font, font.lab, cex, cex.lab, cex.main, lwd, 
              lty1, lty2, pch1, pch2, col1, col2, type="p")

  plotTickLabel(x1,x2,Lx1,Ux1,Lx2,Ux2, axis=1, font.lab, 
                lwd, lty1, lty2, col1, col2) 
  plotTickLabel(tmpy1,tmpy2,Ly1,Uy1,Ly2,Uy2, axis=2, font.lab, 
                lwd, lty1, lty2, col1, col2) 

  invisible(list(sepValx=sepValx, sepValy=sepValy, Q2=Q2))

}

# Obtain first two coordinates of the rotated data sets
# y1, y2 -- clusters 1 and 2
# projDir -- projection direction
# alpha, ITMAX, eps, quiet -- same as other functions
getRotateData<-function(y1, y2, projDir, 
                        iniProjDirMethod=c("SL", "naive"), 
                        projDirMethod=c("newton", "fixedpoint"), 
                        alpha=0.05, ITMAX=10, 
                        eps=1.0e-10, quiet=TRUE)
{ 
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(is.matrix(y1)==TRUE) { n1<-nrow(y1) } 
  else { n1<-1 }
  if(is.matrix(y2)==TRUE) { n2<-nrow(y2) } 
  else { n2<-1 }
  if(n1==1) 
  { mu1<-as.vector(y1)
    p<-length(y1) 
  } 
  else 
  { mu1<-apply(y1, 2, mean, na.rm=TRUE) 
    p<-ncol(y1)
  }
  if(n2==1) { mu2<-as.vector(y2) } 
  else { mu2<-apply(y2, 2, mean, na.rm=TRUE) }
  if(n1==1 || n2==1)
  { cl<-rep(1:2, c(n1, n2))
    y1<-matrix(y1, nrow=n1, ncol=p)
    y2<-matrix(y2, nrow=n2, ncol=p)
    y<-rbind(y1,y2)
    out.label<-NULL
  }
  else  
  { y<-as.matrix(rbind(y1, y2)); cl<-rep(1:2, c(n1, n2)); }
  p<-ncol(y)
  if(p<2)
  { stop("Dimension should be greater than 2!\n") }
  # construct an orthogonal matrix whose first column is projDir 
  Q<-MOrthogonal(projDir)
  # rotate data
  ry<-y%*%Q
  if(p>2)
  { ry1<-ry[,-1] # delete the first dimension
    # we will find a optimal direction in the rest of dimension.
    # the obtained optimal direction will orthogonal with "projDir"
    if(n1>1){tmpy1<-ry1[cl==1,]}else{tmpy1<-matrix(ry1,nrow=1,ncol=p-1)}
    if(n2>1){tmpy2<-ry1[cl==1,]}else{tmpy2<-matrix(ry1,nrow=1,ncol=p-1)}

    tmp<-projDirData(tmpy1, tmpy2, iniProjDirMethod, projDirMethod, 
                alpha, ITMAX, eps, quiet)
    aa<-c(0, tmp$projDir) 

    aa<-as.vector(Q%*%aa) # aa^T a=0
    Q2<-cbind(Q[,1], aa) # data will project into these two dimension
  } 
  else { Q2<-Q }
  x<-y%*%Q2

  return(list(x=x,cl=cl,Q2=Q2))
}


# Visualize clusters by projecting them into a 2-dimensional space.
# The 2-dimensional space is spanned by the first 2 eigenvectors of the
# matrix B:
#\begin{eqnarray*}
# B&=&{2\over
# k_0(k_0-1)}\sum_{i=1}^{k_0-1}\sum_{j=i+1}^{k_0}\E\left[(\Y_i-\Y_j)
# (\Y_i-\Y_j)^T\right]\\
# &=&{2\over k_0}\sum_{i=1}^{k_0}\Sigma_i+{2\over
# k_0(k_0-1)}\sum_{i<j}\left(\thbf_i-\thbf_j\right)
# \left(\thbf_i-\thbf_j\right)^T,
#\end{eqnarray*}
#
# y -- nxp data matrix
# cl -- a partition of the n data points in y
# outlierLabel -- integer (default 0) for detected outliers
# projMethod -- indicates that we use our Eigen method "Eigen" or 
#    method DMS of Dhillon, Modha, Spangler (2002) CSDA v 41, pp 59-90.
# title -- the title of the plot.
viewClusters<-function(y, cl, outlierLabel=0,
  projMethod="Eigen", xlim=NULL, ylim=NULL,
  xlab="1st projection direction", ylab="2nd projection direction", 
  title="Scatter plot of 2-D Projected Clusters",
  font=2, font.lab=2, cex=1.2, cex.lab=1.2) 
{ 
  projMethod<-match.arg(projMethod, choices=c("Eigen", "DMS"))

  y<-as.matrix(y)
  p<-ncol(y);
  y2<-y
  cl2<-cl

  # remove outliers
  out.label<-which(cl==outlierLabel)
  if(length(out.label)>0) 
  { cl2<-cl[-out.label]
    y2<-y[-out.label,,drop=FALSE]
  }

  cl2.u<-unique(cl2)
  k0<-length(cl2.u)

  # project data to the 2-dimensional space spanned by the first 2 
  # eigenvectors of the between cluster distance matrix B.
  if(projMethod=="DMS") { tmp<-DMSProj(y2, cl2, outlierLabel) }
  else { tmp<-eigenProj(y2, cl2, outlierLabel) }
  #ev<-tmp$ev # eigenvalues of the matrix B
  Q<-tmp$Q # the all eigenvectors of B
  B<-tmp$B # the between distance matrix B

  # we only plot the first 2 dimension, using all data points
  # including outliers
  proj<-y%*%Q[,1:2]
 
  if(missing(xlim)) { xlim<-range(proj[,1]) } 
  else { xlim<-xlim }
  if(missing(ylim)) { ylim<-range(proj[,2]) } 
  else { ylim<-ylim }

  # plot projected clusters
  plot(proj[cl==cl2.u[1], 1], proj[cl==cl2.u[1], 2], xlim=xlim, ylim=ylim,
    xlab=xlab,ylab=ylab,col=1, pch=1,font=font,cex=cex,
    font.lab=font.lab,cex.lab=cex.lab)
  # plot outliers
  if(length(out.label)>0)
  { points(proj[cl==outlierLabel,1], proj[cl==outlierLabel,2], col=2, pch=2) }

  title(title,cex.main=1.5); 
  start<-1
  # plot the data points in the remaining clusters
  if(k0>1)
  { for(i in 2:k0)
    { points(proj[cl==cl2.u[i],1], proj[cl==cl2.u[i],2], col=start+i,pch=start+i) }
  }

  res<-list(B=B, Q=Q, proj=proj)
  invisible(res)
}

# Eigen approach to find the projection directions to visualize clusters.
# Project data to the 2-dimensional space spanned by the first two 
# eigenvectors of the between cluster distance matrix B
# B=(2/ k0)*sum_{i=1}^{k_0}Sigmai
#  +(2/(k_0(k_0-1))sum_{i<j}(thetai-thetaj)(thetai-thetaj)^T, 
# where thetai and Sigmai 
# are the mean vector and covariance matrix for i-th cluster.
#
# y -- nxp data points
# cl -- a partition of the n data points in y
# outlierLabel -- integer (default 0) for detected outliers
eigenProj<-function(y,cl, outlierLabel=0)
{ y<-as.matrix(y); p<-ncol(y);
  # remove outliers
  out.label<-which(cl==outlierLabel)
  if(length(out.label)>0) { cl2<-cl[-out.label]; y2<-y[-out.label,]; }
  else { cl2<-cl; y2<-y; }
  cl.set<-sort(unique(cl2))
  k0<-length(cl.set)

  # obtain cluster centers and covariance matrices
  s<-array(0,c(p,p,k0))
  mu.mat<-matrix(0, nrow=k0, ncol=p)
  # W = sum_{i=1}^{k0} Sigma_i /k0
  W<-matrix(0,nrow=p, ncol=p)
  for(i in 1:k0)
  { yi<-y2[cl==cl.set[i],]
    mu.mat[i,]<-apply(yi, 2, mean, na.rm=TRUE)
    s[,,i]<-cov(yi)
    W<-W+s[,,i]
  }
  W<-W/k0

  if(k0>1)
  { # B=sum_{i=1}^{k0-1}sum_{j=(i+1)}^{k0}
    # (theta_i-theta_j)*(theta_i-theta_j)^T / (k0*(k0-1))
    B<-matrix(0, nrow=p, ncol=p)
    for(i in 1:(k0-1))
    { mui<-mu.mat[i,]
      for(j in (i+1):k0)
      { muj<-mu.mat[j,]
        B<-B+outer(mui-muj, mui-muj)
      }
    }
    B<-W+(B/(k0*(k0-1)))
  } else { B<-W }
  B<-2*B

  eg<-eigen(B, symmetric=TRUE)
  Q<-eg$vectors
  # The i-th column of Q is the eigenvector of B
  # corresponding to the i-th eigenvalue
  res<-list(Q=Q, B=B)
  return(res)
}

# Dhillon et al.'s (2002) approach to find the projection directions to 
# visualize clusters.
#  Dhillon I. S., Modha, D. S. and Spangler, W. S. (2002)
#  Class visualization of high-dimensional data with applications.
#  \emph{computational Statistics and Data Analysis}, \bold{41}, 59--90.
# 
# project data to the t-dimensional space spanned by the first two eigen
# vectors of the between cluster distance matrix B
# B=\sum_{i=2}^{k_0}\sum_{j=1}^{i-1}
#   n_in_j(\theta_i-\theta_j)(\theta_i-\theta_j)^T
#
# y -- nxp data points
# cl -- a partition of the n data points in y
# outlierLabel -- integer (default 0) for detected outliers
DMSProj<-function(y, cl, outlierLabel=0)
{ y<-as.matrix(y)
  p<-ncol(y);
  # remove outliers
  out.label<-which(cl==outlierLabel)
  if(length(out.label)>0) { cl2<-cl[-out.label]; y2<-y[-out.label,]; }
  else { cl2<-cl; y2<-y; }
  cl.set<-sort(unique(cl2))
  k0<-length(cl.set)
  n.set<-tapply(rep(1,length(cl2)), cl2, sum, na.rm=TRUE)

  # obtain cluster centers and covariance matrices
  mu.mat<-matrix(0, nrow=k0, ncol=p)
  for(i in 1:k0)
  { yi<-y2[cl==cl.set[i],]
    if(n.set[i]>1)
    { mu.mat[i,]<-apply(yi, 2, mean, na.rm=TRUE) }
    else { mu.mat[i,]<-yi }
  }

  if(k0>1)
  { # B=\sum_{i=2}^{k_0}\sum_{j=1}^{i-1}
    # n_in_j(\theta_i-\theta_j)(\theta_i-\theta_j)^T
    B<-matrix(0, nrow=p, ncol=p)
    for(i in 2:k0)
    { mui<-as.vector(mu.mat[i,])
      for(j in 1:(i-1))
      { muj<-as.vector(mu.mat[j,])
        ni<-n.set[i]; nj<-n.set[j];
        tem1<-(mui-muj)
        tem2<-tem1%*%t(tem1)
        B<-B+n.set[i]*n.set[j]*tem2
      }
    }
  }
  else { stop("the number of clusters is 0!\n") }
  eg<-eigen(B, symmetric=TRUE)
  Q<-eg$vectors
  # The i-th column of Q is the eigenvector of B
  # corresponding to the i-th eigenvalue
  res<-list(Q=Q, B=B)
  return(res)
}


