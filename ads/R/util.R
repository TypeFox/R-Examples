overlapping.polygons<-function(listpoly) {
	stopifnot(unlist(lapply(listpoly,is.poly)))
	nbpoly<-length(listpoly)
	res<-rep(FALSE,(nbpoly-1))
	for(i in 1:(nbpoly-1)) {
		for(j in (i+1):nbpoly) {
		 if(abs(overlap.poly(listpoly[[i]],listpoly[[j]]))>.Machine$double.eps^0.5)
			res[i]<-TRUE
		}
	}
	return(res)
}

#from area.xypolygon{spatstat}
#return area>0 when (xp,yp) vertices are ranked anticlockwise
#or area<0 when (xp,yp) vertices are ranked clockwise
area.poly<-function(xp,yp) {
    nedges <- length(xp)
    yp <- yp - min(yp)
    nxt <- c(2:nedges, 1)
    dx <- xp[nxt] - xp
    ym <- (yp + yp[nxt])/2
    -sum(dx * ym)
}

in.circle<-function(x,y,x0,y0,r0,bdry=TRUE) {
	stopifnot(length(x)==length(y))
	l<-length(x)
	inside<-vector(mode="logical",length=l)
	for(i in 1:l) {
		if(bdry) {
			if(((x[i]-x0)^2+(y[i]-y0)^2)<=(r0^2))
				inside[i]<-TRUE
		}
		else {
			if(((x[i]-x0)^2+(y[i]-y0)^2)<(r0^2))
				inside[i]<-TRUE
		}
	}
	return(inside)
}

in.rectangle<-function(x,y,xmin,ymin,xmax,ymax,bdry=TRUE) {
	stopifnot(length(x)==length(y))
	stopifnot((xmax-xmin)>0)
	stopifnot((ymax-ymin)>0)
	rect<-list(x=c(xmin,xmax,xmax,xmin),y=c(ymin,ymin,ymax,ymax))
	return(in.poly(x,y,rect,bdry))
}

in.triangle<-function(x,y,ax,ay,bx,by,cx,cy,bdry=TRUE) {
	stopifnot(length(x)==length(y))
	tri<-list(x=c(ax,bx,cx),y=c(ay,by,cy))
	return(in.poly(x,y,tri,bdry))
}

#modified from plot.ppp{spatstat}
adjust.marks.size<-function(marks,window,maxsize=NULL) {
	if(is.null(maxsize)) {
		if("rectangle"%in%window$type)
			diam<-sqrt((window$xmin-window$xmax)^2+(window$ymin-window$ymax)^2)
		else if("circle"%in%window$type)
			diam<-2*window$r0
		maxsize<-min(1.4/sqrt(pi*length(marks)/area.swin(window)),diam*0.07)
	}
	mr<-range(c(0,marks))
	maxabs<-max(abs(mr))
	if(diff(mr)<4*.Machine$double.eps||maxabs<4*.Machine$double.eps) {
		ms<-rep(0.5*maxsize,length(marks))
		mp.value<-mr[1]
		mp.plotted<-0.5*maxsize
	}
	else {
		scal<-maxsize/maxabs
		ms<-marks*scal
		mp.value<-pretty(mr)
		mp.plotted<-mp.value*scal
	}
	return(ms)
}

#modified from verify.xypolygon{spatstat}
is.poly<-function (p) {
	stopifnot(is.list(p))
	stopifnot(length(p)==2,length(names(p))==2)
	stopifnot(identical(sort(names(p)),c("x","y")))
	stopifnot(!is.null(p$x),!is.null(p$y))
	stopifnot(is.numeric(p$x),is.numeric(p$y))
	stopifnot(length(p$x)==length(p$y))
	return(TRUE)
}

testInteger<-function(i) {
	if(as.integer(i)!=i) {
		warning(paste(substitute(i),"=",i," has been converted to integer : ",as.integer(i),sep=""),call.=FALSE)
		i<-as.integer(i)
	}
	return(i)
}

testIC<-function(nbSimu,lev) {
	if(lev*(nbSimu+1)<5) {
		warning(paste(
			"Low validity test: a*(n+1) < 5\n  Significance level: a = ",lev,
			"\n  Number of simulations: n = ",nbSimu,"\n",sep=""))
	}
}

#from spatstat
overlap.poly <- function(P, Q) {
  # compute area of overlap of two simple closed polygons 
 # verify.xypolygon(P)
  #verify.xypolygon(Q)
  
  xp <- P$x
  yp <- P$y
  np <- length(xp)
  nextp <- c(2:np, 1)

  xq <- Q$x
  yq <- Q$y
  nq <- length(xq)
  nextq <- c(2:nq, 1)

  # adjust y coordinates so all are nonnegative
  ylow <- min(c(yp,yq))
  yp <- yp - ylow
  yq <- yq - ylow

  area <- 0
  for(i in 1:np) {
    ii <- c(i, nextp[i])
    xpii <- xp[ii]
    ypii <- yp[ii]
    for(j in 1:nq) {
      jj <- c(j, nextq[j])
      area <- area +
        overlap.trapez(xpii, ypii, xq[jj], yq[jj])
    }
  }
  return(area)
}

#from spatstat
overlap.trapez <- function(xa, ya, xb, yb, verb=FALSE) {
  # compute area of overlap of two trapezia
  # which have same baseline y = 0
  #
  # first trapezium has vertices
  # (xa[1], 0), (xa[1], ya[1]), (xa[2], ya[2]), (xa[2], 0).
  # Similarly for second trapezium
  
  # Test for vertical edges
  dxa <- diff(xa)
  dxb <- diff(xb)
  if(dxa == 0 || dxb == 0)
    return(0)

  # Order x coordinates, x0 < x1
  if(dxa > 0) {
    signa <- 1
    lefta <- 1
    righta <- 2
    if(verb) cat("A is positive\n")
  } else {
    signa <- -1
    lefta <- 2
    righta <- 1
    if(verb) cat("A is negative\n")
  }
  if(dxb > 0) {
    signb <- 1
    leftb <- 1
    rightb <- 2
    if(verb) cat("B is positive\n")
  } else {
    signb <- -1
    leftb <- 2
    rightb <- 1
    if(verb) cat("B is negative\n")
  }
  signfactor <- signa * signb # actually (-signa) * (-signb)
  if(verb) cat(paste("sign factor =", signfactor, "\n"))

  # Intersect x ranges
  x0 <- max(xa[lefta], xb[leftb])
  x1 <- min(xa[righta], xb[rightb])
  if(x0 >= x1)
    return(0)
  if(verb) {
    cat(paste("Intersection of x ranges: [", x0, ",", x1, "]\n"))
    abline(v=x0, lty=3)
    abline(v=x1, lty=3)
  }

  # Compute associated y coordinates
  slopea <- diff(ya)/diff(xa)
  y0a <- ya[lefta] + slopea * (x0-xa[lefta])
  y1a <- ya[lefta] + slopea * (x1-xa[lefta])
  slopeb <- diff(yb)/diff(xb)
  y0b <- yb[leftb] + slopeb * (x0-xb[leftb])
  y1b <- yb[leftb] + slopeb * (x1-xb[leftb])
  
  # Determine whether upper edges intersect
  # if not, intersection is a single trapezium
  # if so, intersection is a union of two trapezia

  yd0 <- y0b - y0a
  yd1 <- y1b - y1a
  if(yd0 * yd1 >= 0) {
    # edges do not intersect
    areaT <- (x1 - x0) * (min(y1a,y1b) + min(y0a,y0b))/2
    if(verb) cat(paste("Edges do not intersect\n"))
  } else {
    # edges do intersect
    # find intersection
    xint <- x0 + (x1-x0) * abs(yd0/(yd1 - yd0))
    yint <- y0a + slopea * (xint - x0)
    if(verb) {
      cat(paste("Edges intersect at (", xint, ",", yint, ")\n"))
      points(xint, yint, cex=2, pch="O")
    }
    # evaluate left trapezium
    left <- (xint - x0) * (min(y0a, y0b) + yint)/2
    # evaluate right trapezium
    right <- (x1 - xint) * (min(y1a, y1b) + yint)/2
    areaT <- left + right
    if(verb)
      cat(paste("Left area = ", left, ", right=", right, "\n"))    
  }

  # return area of intersection multiplied by signs 
  return(signfactor * areaT)
}

#TRUE: les points sur la bordure sont = inside
in.poly<-function(x,y,poly,bdry=TRUE) {
	stopifnot(is.poly(poly))
	xp <- poly$x
	yp <- poly$y
	npts <- length(x)
	nedges <- length(xp)   # sic

  score <- rep(0, npts)
  on.boundary <- rep(FALSE, npts)
	temp <- .Fortran(
		"inpoly",
		x=as.double(x),
		y=as.double(y),
		xp=as.double(xp),
		yp=as.double(yp),
		npts=as.integer(npts),
		nedges=as.integer(nedges),
		score=as.double(score),
		onbndry=as.logical(on.boundary),
		PACKAGE="ads"
	)
	score <- temp$score
	on.boundary <- temp$onbndry
	score[on.boundary] <- 1
	res<-rep(FALSE,npts)
	res[score==(-1)]<-TRUE
	if(bdry)
		res[score==1]<-TRUE
	return(res)
}

####################
convert<-function(x) {
r<-alist()
x<-as.matrix(x)
for(i in 1:dim(x)[1])
	r[[i]]<-data.frame(x=c(x[i,1],x[i,5],x[i,3]),y=c(x[i,2],x[i,6],x[i,4]))
return(r)
}

convert2<-function(x){ # liste de liste vers df 6 var
x<-unlist(x)
mat<-matrix(x,ncol=6,byrow=TRUE,dimnames=list(NULL,c("ax","bx","cx","ay","by","cy")))
#mat<-cbind(tmp$ax,tmp$ay,tmp$bx,tmp$by,tmp$cx,tmp$cy)
r<-data.frame(ax=mat[,1],ay=mat[,4],bx=mat[,2],by=mat[,5],cx=mat[,3],cy=mat[,6])
return(r)
}


read.tri<-function(X) {
	res<-NULL
	tabtri<-read.table(X)	
	
	n<-length(tabtri[,1])
	
	for(i in 1:n) {
		tri<-list(x=c(tabtri[,1][i],tabtri[,3][i],tabtri[,5][i]),
				  y=c(tabtri[,2][i],tabtri[,4][i],tabtri[,6][i]))
				  
		#if(area.xypolygon(tri)>0){
		if(area.poly(tri$x,tri$y)>0){
			tri<-list(x=c(tabtri[,5][i],tabtri[,3][i],tabtri[,1][i]),
				  y=c(tabtri[,6][i],tabtri[,4][i],tabtri[,2][i]))
		}
				  
		res<-c(res,list(tri))
	}
	if(length(res)==1) {
		res<-unlist(res,recursive=FALSE)
	}
	return(res)
}

transpose<-function(x,y) {

	nbTri<-length(x)/3
	
	res<-.C("transpose",x=as.double(x),y=as.double(y),nbTri=as.integer(nbTri),
			x1=double(nbTri),y1=double(nbTri),x2=double(nbTri),y2=double(nbTri),
			x3=double(nbTri),y3=double(nbTri),PACKAGE="ads")
		
	list(x1=res$x1,y1=res$y1,x2=res$x2,y2=res$y2,x3=res$x3,y3=res$y3)
}
##############
#subsetting dist objects
#sub is a logical vector of True/False
subsetdist<-function(dis,sub) {
	mat<-as.matrix(dis)
	k<-dimnames(mat)[[1]]%in%sub
	submat<-mat[k,k]
	return(as.dist(submat))
}
	
#ordering dist objetcs on ind
sortmat<-function(dis,ind) {
	mat<-as.matrix(dis)[,ind]
	mat<-mat[ind,]
	return(as.dist(mat))
}


