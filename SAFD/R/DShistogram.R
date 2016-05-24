DShistogram <-
function(XX,limx=NA,npart=10,nl=101,pic=TRUE,pdf=FALSE){
 #function makes a partition of the interval stated in xlim of nbreaks elements
 #XX...fuzzy sample (list as always)
 #xlim...limits of the histrogram - if NA then the max and min of the supps will be taken
 #nl...number of levels
 #make use of frequency function
 #construct 2 dim matrix and use 3d plot via persp and contour plot
 
 #construct limits
 if(length(limx)<=1|limx[2]<=limx[1]){
  a<-XX[[1]]$x[1]
  b<-XX[[1]]$x[nrow(XX[[1]])]
  if(length(XX)>1){
    for (i in 2:length(XX)){
     a<-min(a,XX[[i]]$x[1])
     b<-max(b,XX[[i]]$x[nrow(XX[[i]])])
    }
   }
  limx<-c(a,b) 
 }
 
 k<-length(XX)
  if(k>500){
    ygrid<-seq(0,1,length=501)
  }
  if(k<=500){
  ygrid<-sort(union(seq(0,1,length=(k+1)),seq(0,1,length=101)))
  }
 
 breaks<-seq(limx[1],limx[2],length=npart+1)
 FR<-vector("list",length=npart)
 FR2<-vector("list",length=npart)
 for (i in 1:npart){
  FR[[i]]<-DSfrequency(XX,breaks[i:(i+1)],0,nl)
  print(i)
  R<-FR[[i]][(nl+1):(2*nl),]
  a<-approx(R$x,R$alpha,xout=ygrid,yleft=R$alpha[1],yright=R$alpha[nl],
   method="constant",f=1,ties="max")
  L<-FR[[i]][1:nl,]
  b<-approx(L$x,L$alpha,xout=ygrid,yleft=L$alpha[1],yright=L$alpha[nl],
   method="constant",f=0,ties="max") 
  
  value<-ifelse(a$y>=b$y,b$y,a$y) 
  FR2[[i]]<-data.frame(x=ygrid,y=value)
  }
 #construct grid for y-coordinate in plotting
 grid1<-breaks+(breaks[2]-breaks[1])/1000
 grid2<-breaks-(breaks[2]-breaks[1])/1000
 grid3<-c(grid1,grid2)
 grid3<-sort(subset(grid3,grid3>=min(breaks)&grid3<=max(breaks)))
 gridx<-grid3
 
 gridy<-ygrid
 M<-matrix(numeric(npart*length(gridy)),ncol=length(gridy))
 for (i in 1:npart){
  M[i,]<-FR2[[i]]$y
 }
 
 M2<-M[rep(1:npart, rep(2,npart)),]

  k<-length(XX)
  lower<-rep(0,k)
  upper<-lower
  for (j in 1:k){
   lower[j]<-min(XX[[j]]$x)
   upper[j]<-max(XX[[j]]$x)
   }
  lim_temp<-c(min(lower),max(upper))
  
 if(pdf==TRUE){
  pdf(file="histo.pdf",width=12,height=8)
  #BBreaks<-list(length=length(breaks))
  #for (m in 1:length(breaks)){
  # BBreaks[[m]]<-data.frame(x=rep(breaks[m],2),alpha=c(-0.05,1.05))
  #}
  #plot(XX[[1]],type="l", xlim=lim_temp,lwd=0.3,xlab=" ", ylab=" ",cex.main=1, col="gray50",
  #  main=paste("Sample",sep=""))
  #for (j in 2:min(k,200)){
  # lines(XX[[j]],type="l",lwd=0.3,col="gray50")
  #}
  #for (m in 1:length(breaks)){
  # lines(BBreaks[[m]],type="l",col="red",lwd=2)
  # }
   
 color<-rainbow(100,start=.7,end=.17)
  # Compute the z-value at the facet centres
 zfacet <- M2[-1, -1] + M2[-1, -ncol(M2)] + M2[-nrow(M2), -1] + M2[-nrow(M2), -ncol(M2)]
 facetcol <- cut(zfacet, 100)
 M<-M2
 
 #calculate plot limit for y-coordinate
 colmax<-rep(0,trunc(length(gridy)/10))
  for (i in 1:trunc(length(gridy)/10)){
 colmax[i]<-max(M[,10*i])
 }
 Cut<-data.frame(nr=seq(1,length(colmax),by=1),colmax=colmax)
 Cut<-subset(Cut,Cut$colmax>0)
 cutindex<-min(round(10*Cut$nr[nrow(Cut)]*1.25,0),length(gridy))
 ym<-min(gridy[10*Cut$nr[nrow(Cut)]]*1.25,1)
 #print(ym)
 
 Mp<-M[,1:cutindex]
 gridyp<-gridy[1:cutindex]
    
 persp(gridx,gridyp,Mp, xlab="x",  ylab="upper/lower frequency", zlab=expression(alpha),
    xlim=limx, main=paste("Histogram 3d",sep=""),cex.main=1,
    theta = -45, phi = 35, expand = 0.35, col=color[facetcol],
    shade = 0.25, ticktype = "detailed",border=NA)
 persp(gridx,gridyp,Mp, xlab="x",  ylab="upper/lower frequency", zlab=expression(alpha),
    xlim=limx, main=paste("Histogram 3d",sep=""),cex.main=1,
    theta = 45, phi = 35, expand = 0.35, col=color[facetcol],
    shade = 0.25, ticktype = "detailed",border=NA)

 image(gridx,gridyp,Mp, xlab="x",  ylab="upper/lower frequency", xlim=limx,
          col=rainbow(100,start=.7,end=.17),cex.axis=1,
          main=paste("Histogram level view","\n",
            "(black lines denote 1-cut, white lines 0.5-cut)",sep=""),cex.main=1)
 contour(gridx,gridyp,Mp, xlab=NA,  ylab=NA, xlim=limx,lwd=c(1.5,1.5),
        levels = seq(0.5,1,by=0.5), add = TRUE, col = c("white","black"),
         lty = c(1,1), drawlabels=FALSE)

 dev.off()
 }
 if(pic==TRUE){
   color<-rainbow(100,start=.7,end=.17)
   # Compute the z-value at the facet centres
  zfacet <- M2[-1, -1] + M2[-1, -ncol(M2)] + M2[-nrow(M2), -1] + M2[-nrow(M2), -ncol(M2)]
  facetcol <- cut(zfacet, 100)
  M<-M2
  
  #calculate plot limit for y-coordinate
  colmax<-rep(0,trunc(length(gridy)/10))
   for (i in 1:trunc(length(gridy)/10)){
  colmax[i]<-max(M[,10*i])
  }
  Cut<-data.frame(nr=seq(1,length(colmax),by=1),colmax=colmax)
  Cut<-subset(Cut,Cut$colmax>0)
  cutindex<-min(round(10*Cut$nr[nrow(Cut)]*1.25,0),length(gridy))
  ym<-min(gridy[10*Cut$nr[nrow(Cut)]]*1.25,1)
  #print(ym)
  
  Mp<-M[,1:cutindex]
  gridyp<-gridy[1:cutindex]
      
  persp(gridx,gridyp,Mp, xlab="x",  ylab="upper/lower frequency", zlab=expression(alpha),
    xlim=limx, main=paste("Histogram 3d",sep=""),cex.main=1,
    theta = -45, phi = 35, expand = 0.35, col=color[facetcol],
    shade = 0.25, ticktype = "detailed",border=NA)
  dev.new()  
  image(gridx,gridyp,Mp, xlab="x",  ylab="upper/lower frequency", xlim=limx,
          col=rainbow(100,start=.7,end=.17),cex.axis=1,
          main=paste("Histogram level view","\n",
            "(black lines denote 1-cut, white lines 0.5-cut)",sep=""),cex.main=1)
 contour(gridx,gridyp,Mp, xlab="",  ylab="", xlim=limx,lwd=c(1.5,1.5),
        levels = seq(0.5,1,by=0.5), add = TRUE, col = c("white","black"),
         lty = c(1,1), drawlabels=FALSE)

 }
 H<-list(gridx=gridx,gridy=gridy,M=M,breaks=breaks)
invisible(H)
}
