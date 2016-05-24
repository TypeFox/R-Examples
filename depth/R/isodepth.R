isodepth=function(x,dpth=NULL,output=FALSE,twodim=TRUE,mustdith=FALSE,maxdith=50,dithfactor=10,trace.errors=TRUE,eps=1e-8,factor=.8,xlab="X",ylab="Y",zlab="Tukey's depth",colcontours=NULL,...){

  if(is.data.frame(x)) x=as.matrix(x)
  if(is.list(x)) {
    m=length(x)
    n=length(x[[1]])
    y=matrix(0,n,m)
    for(i in 1:m){
      y[,i]=x[[i]]
      if(length(x[[i]])!=n){ stop("When using a list, each element must be a vector of the same length.") }
    }
    x=y
  }

  p=length(x[1,])
  n=length(x[,1])

  if(p>n) { warning(message=paste("Is your data ",n," points in ",p," dimensions.\nIf not, you should transpose your data matrix.")) }
  if(p!=2) { stop("Data must be bivariate.\n") }
  
  	if(is.null(dpth)) { dpth=seq(1,floor(n/2)) }

	ndpth=length(dpth)
	y=x[,2]
	x=x[,1]

	maxnum=floor(4*n*sqrt(n)+1)

	zz=.Fortran("halfmed",
			as.numeric(x),
			as.numeric(y),
			as.integer(n),
			integer(1),
			numeric(2),
			xcont=numeric(n*(n-1)/2),
			ycont=numeric(n*(n-1)/2),
			ncont=integer(ndpth),
			as.integer(dpth),
			as.integer(ndpth),
			as.integer(1),
			as.integer(maxnum),
			err=integer(1),
			as.numeric(eps),
			as.numeric(dithfactor),
			as.integer(maxdith),
			as.integer(mustdith),
			missing=integer(ndpth),
			as.numeric(factor),
			PACKAGE="depth")

	if(zz$err==-1&&mustdith==FALSE){ stop("Points are not in general position.") }

	missed=-1
	needdith=-1
	toomuchdith=-1

	if(zz$missing[1]==-1){missed=dpth[1]}
	if(zz$missing[1]==-2){needdith=dpth[1]}
	if(zz$missing[1]==-3){toomuchdith=dpth[1]}
	if(ndpth>1){
		for(i in 2:ndpth){
			if(zz$missing[i]==-1){
				if(missed[1]==-1){
					missed=dpth[i]}
				else { missed=c(missed,dpth[i])}
			}
			if(zz$missing[i]==-2){
				if(needdith[1]==-1){
					needdith=dpth[i] }
				else { needdith=c(needdith,dpth[i])}
			}
			if(zz$missing[i]==-3){
				if(toomuchdith[1]==-1) {
					toomuchdith=dpth[i] }
				else { toomuchdith=c(toomuchdith,dpth[i])}
			}
		}
	}
	if(missed[1]!=-1) { warning("Depth contours ",paste(missed,collapse=",")," do not exist.") }
	if(needdith[1]!=-1) { warning("Numerical problems encountered for contours ",needdith,". They may be wrong.\n") }
	if(toomuchdith[1]!=-1) { cat("\nDepths ",toomuchdith," could not be calculated even after ",maxdith," steps of ventilation.\n") }

	flag=1
	for(i in 1:ndpth) {
		t=zz$ncont[i]
		for(j in 1:t) {
			if(is.na(zz$xcont[flag])||is.na(zz$ycont[flag])){
				zz$ncont[i]=zz$ncont[i]-1
				zz$xcont=c(zz$xcont[1:flag-1],zz$xcont[flag+1:length(zz$xcont)])
				zz$ycont=c(zz$ycont[1:flag-1],zz$ycont[flag+1:length(zz$ycont)])
			flag=flag-1
			}

			flag=flag+1
		}
	}

  if(twodim==TRUE & output==FALSE){
    
        if(is.null(col)){  col="black" }
        if(is.null(colcontours)){  colcontours=rep("black",ndpth) }
	if(length(colcontours)<ndpth){ colcontours[1:ndpth]=colcontours }
		
      	plot(x,y,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),xlab=xlab,ylab=ylab,...)
	nc=0
	for (i in 1:ndpth) {
		xx=0
		yy=0
		if (zz$ncont[i]>=1) {
			xx=zz$xcont[nc+1]
			yy=zz$ycont[nc+1]
		
			for (j in 2:zz$ncont[i]) {
				xx=c(xx,zz$xcont[nc+j])
				yy=c(yy,zz$ycont[nc+j])
			}
		}
		if (zz$ncont[i]>2) {
			xxx=c(zz$xcont[nc+1],zz$xcont[nc+zz$ncont[i]])
			yyy=c(zz$ycont[nc+1],zz$ycont[nc+zz$ncont[i]])

		if(trace.errors == TRUE || zz$missing[i] >= 0) {
			lines(xxx,yyy,col=colcontours[i])
		}

		}
		nc=nc+zz$ncont[i]

		if (zz$ncont[i]>2){
			if(trace.errors == TRUE || zz$missing[i] >=0 ) {	
		  lines(xx,yy,col=colcontours[i])
		}
		}
		if(zz$ncont[i]==1) {
			if(trace.errors == TRUE || zz$missing[i] >= 0) {	
				points(xx,yy,pch=3,col=colcontours[i])
			}
		}
		
	}
  
    }
    
    if(twodim==FALSE & output==FALSE){
 
 	require(rgl)
	xx=NULL
	yy=NULL
	Dpth=NULL
	ncont=NULL
	nc=0
	for (i in 1:ndpth) {
		if (zz$ncont[i]>=1) {
			if(trace.errors == TRUE || zz$missing[i] >= 0) {	
				xx=c(xx,zz$xcont[nc+1])
				yy=c(yy,zz$ycont[nc+1])
			
				if (zz$ncont[i]>2) {

					for (j in 2:zz$ncont[i]) {
						xx=c(xx,zz$xcont[nc+j])
						yy=c(yy,zz$ycont[nc+j])
					}
				}

				ncont=c(ncont,zz$ncont[i])
				Dpth=c(Dpth,dpth[i])

			}
			nc=nc+zz$ncont[i]
		}

	}
  	plot3d(x,y,rep(0,length(x)),type="p",size=3,zlim=c(0,max(Dpth)),xlab=xlab,ylab=ylab,zlab=zlab,...)
	  ncont2=c(0,cumsum(ncont))
	  Dpth=c(0,Dpth)
	  for(i in 1:length(ncont)){
	    keep=((ncont2[i])+1):(ncont2[i+1])
	    px=c(rep(xx[keep],each=2),xx[keep[1]],xx[keep[1]])
	    py=c(rep(yy[keep],each=2),yy[keep[1]],yy[keep[1]])
	    pz1=rep(Dpth[c(i,(i+1),i+1,i)],length(keep))[1:length(px)] 
	    pz2=rep(Dpth[c(i+1,(i),i,i+1)],length(keep))[1:length(px)]
	    plot3d(px,py,pz1,type='l',add=TRUE)
	    plot3d(px,py,pz2,type='l',add=TRUE)    
	  }
    }
    
    if(output==TRUE){
    	rep=vector("list",ndpth)
 	nc=0
	for (i in 1:ndpth) {
		xx=0
		yy=0
		if (zz$ncont[i]>=1) {
			rep[[i]]=cbind(zz$xcont[(nc+1):(nc+zz$ncont[i])],zz$ycont[(nc+1):(nc+zz$ncont[i])])
		}
		nc=nc+zz$ncont[i]
	}
      names(rep)=paste("Contour",dpth,sep="")
      return(rep)
    }   
}
