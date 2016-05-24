responseMap <- function(model,data,irreg=FALSE,interpolate=TRUE,margins=TRUE,logscale=TRUE,zlims=NULL,raster=TRUE) UseMethod("responseMap")
responseMap.default <- function(model,data,irreg=FALSE,interpolate=TRUE,margins=TRUE,logscale=TRUE,zlims=NULL,raster=TRUE) {
	concs <- model
	act <- data
	if (ncol(concs)<2) { stop("Argument 'concs' must have at least two columns.") }
	if (is.null(act) & ncol(concs)<3) { stop("If argument 'act' is null, argument 'concs' must have at least 3 columns.") }
	conc1 <- as.vector(concs[,1])
	conc2 <- as.vector(concs[,2])
	if (is.null(act)) { act <- as.vector(concs[,3]) }
	else {
		act <- as.vector(act)
		if (length(act)!=length(conc1)) { stop("Length of argument 'act' must be equal to number of rows in argument 'conc'.") }
	}
	return(makeResponseMap(conc1,conc2,act,irreg,interpolate,margins,logscale,zlims,raster))
}
responseMap.formula <- function(model,data,...) {
	mf <- model.frame(formula=model, data=data)
	concs <- model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(concs,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { concs <- concs[,-i] }
	}
	act <- model.response(mf)
	return(responseMap.default(concs,act,...))
}
makeResponseMap <- function(conc1,conc2,act,irreg=FALSE,interpolate=TRUE,margins=TRUE,logscale=TRUE,zlims=NULL,raster=TRUE) {
	if (!is.null(zlims)) {
		if (length(zlims)!=2) { stop("Argument 'zlims' must be a vector of length 2.") }
		act <- pmin(pmax(act,zlims[1]),zlims[2])
	}
	
	if (logscale) { allpos <- conc1>0 & conc2>0 }
	else { 
		allpos <- rep(TRUE,length(conc1))
		margins <- FALSE
	}
	mainpts <- cbind(conc1[allpos],conc2[allpos])
	mainact <- act[allpos]
	if (logscale) { 
		mainpts <- log(mainpts)
		if (margins) {
			margxpts <- conc1[conc1>0 & conc2==0]
			if (length(margxpts)>0) { 
				margxact <- act[conc1>0 & conc2==0]
				margxpts <- log(margxpts)
			}
			margypts <- conc2[conc1==0 & conc2>0]
			if (length(margypts)>0) {
				margyact <- act[conc1==0 & conc2>0]
				margypts <- log(margypts)
			}
		}
	}
	if (interpolate) {
		nval <- 100
		ix <- seq(min(mainpts[,1]),max(mainpts[,1]),length=nval)
		iy <- seq(min(mainpts[,2]),max(mainpts[,2]),length=nval)
		px <- rep(ix,each=length(iy))
		py <- rep(iy,times=length(ix))
		margxwidth <- margywidth <- 0
		if (margins) { 
			mainywidth <- (max(iy)-min(iy))*nval/(nval-1)
			mainxwidth <- (max(ix)-min(ix))*nval/(nval-1)
			margxwidth <- mainywidth/14
			margywidth <- mainxwidth/14
			margbarwidth <- max(margxwidth,margywidth)
			margxposy <- min(iy)-mainywidth/(2*nval)-margbarwidth
			margyposx <- min(ix)-mainxwidth/(2*nval)-margbarwidth
			margbarwidth <- log10(exp(margbarwidth))
		}
		if (irreg) {
			pz <- linearInterpIrreg2d(mainpts[,1],mainpts[,2],mainact,px,py)
			if (margins && length(margxpts)>0) { margxdf <- data.frame(x=ix,y=margxposy,z=linearInterpIrreg1d(margxpts,margxact,ix)) }
			if (margins && length(margypts)>0) { margydf <- data.frame(x=margyposx,y=iy,z=linearInterpIrreg1d(margypts,margyact,iy)) }
		} else {
			pz <- gaussInterp2d(mainpts[,1],mainpts[,2],mainact,px,py)
			if (margins && length(margxpts)>0) { margxdf <- data.frame(x=ix,y=margxposy,z=gaussInterp1d(margxpts,margxact,ix)) }
			if (margins && length(margypts)>0) { margydf <- data.frame(x=margyposx,y=iy,z=gaussInterp1d(margypts,margyact,iy)) }
		}
		maindf <- data.frame(x=px,y=py,z=pz)
	} else {
		if (irreg) {
			umainpts <- unique(mainpts)
			uz <- rep(0,nrow(umainpts))
			for (i in 1:nrow(umainpts)) { uz[i] <- mean(mainact[mainpts[,1]==umainpts[i,1] & mainpts[,2]==umainpts[i,2]]) }
			maindf <- voronoi(umainpts)
			maindf$z <- uz[maindf$poly]
			if (margins) {
				mainywidth <- max(maindf$y)-min(maindf$y)
				mainxwidth <- max(maindf$x)-min(maindf$x)
				minyval <- min(maindf$y)
				minxval <- min(maindf$x)
				maxyval <- max(maindf$y)
				maxxval <- max(maindf$x)
			}
		} else {
			uxvals <- sort(unique(mainpts[,1]))
			uyvals <- sort(unique(mainpts[,2]))
			uz <- rep(0,length(uxvals)*length(uyvals))
			for (i in 1:length(uxvals)) {
				for (j in 1:length(uyvals)) {
					rel <- mainpts[,1]==uxvals[i] & mainpts[,2]==uyvals[j]
					if (any(rel)) { uz[(i-1)*length(uyvals)+j] <- mean(mainact[rel]) }
					else {
						if (i==1) { relx <- mainpts[,1]==uxvals[i+1] }
						else if (i==length(uxvals)) { relx <- mainpts[,1]==uxvals[i-1] }
						else { relx <- mainpts[,1]==uxvals[i+1] | mainpts[,1]==uxvals[i-1] }
						if (j==1) { rely <- mainpts[,2]==uyvals[j+1] }
						else if (j==length(uyvals)) { rely <- mainpts[,2]==uyvals[j-1] }
						else { rely <- mainpts[,2]==uyvals[j+1] | mainpts[,2]==uyvals[j-1] }
						rel <- (mainpts[,1]==uxvals[i] & rely) | (relx & mainpts[,2]==uyvals[j])
						if (any(rel)) { uz[(i-1)*length(uyvals)+j] <- mean(mainact[rel]) }
						else { uz[(i-1)*length(uyvals)+j] <- NA }
					}
				}
			}
			maindf <- data.frame(x=rep(uxvals,each=length(uyvals)),y=rep(uyvals,times=length(uxvals)),z=uz)
			if (margins) {
				mainywidth <- (max(uyvals)-min(uyvals))*length(uyvals)/(length(uyvals)-1)
				mainxwidth <- (max(uxvals)-min(uxvals))*length(uxvals)/(length(uxvals)-1)
				minyval <- min(uyvals)-mainywidth/(2*length(uyvals))
				minxval <- min(uxvals)-mainxwidth/(2*length(uxvals))
				maxyval <- max(uyvals)+mainywidth/(2*length(uyvals))
				maxxval <- max(uxvals)+mainxwidth/(2*length(uxvals))
			}
		}
		if (margins) {
			margxwidth <- mainywidth/14
			margywidth <- mainxwidth/14
			margbarwidth <- max(margxwidth,margywidth)
			margxposy <- minyval-margbarwidth
			margyposx <- minxval-margbarwidth
			margbarwidth <- log10(exp(margbarwidth))
			if (length(margxpts>0)) {
				umxv <- unique(margxpts)
				umxv <- sort(umxv[umxv>minxval&umxv<maxxval])
				bnds <- (umxv[1:(length(umxv)-1)]+umxv[2:length(umxv)])/2
				bnds <- c(minxval,bnds,maxxval)
				px <- (bnds[1:(length(bnds)-1)]+bnds[2:length(bnds)])/2
				pw <- bnds[2:length(bnds)]-bnds[1:(length(bnds)-1)]
				pw <- log10(exp(pw))
				pz <- rep(0,length(px))
				for (i in 1:length(px)) { pz[i] <- mean(margxact[margxpts==umxv[i]]) }
				margxdf <- data.frame(x=px,y=margxposy,z=pz,w=pw)
			}
			if (length(margypts>0)) {
				umyv <- unique(margypts)
				umyv <- sort(umyv[umyv>minyval&umyv<maxyval])
				bnds <- (umyv[1:(length(umyv)-1)]+umyv[2:length(umyv)])/2
				bnds <- c(minyval,bnds,maxyval)
				py <- (bnds[1:(length(bnds)-1)]+bnds[2:length(bnds)])/2
				pw <- bnds[2:length(bnds)]-bnds[1:(length(bnds)-1)]
				pw <- log10(exp(pw))
				pz <- rep(0,length(py))
				for (i in 1:length(py)) { pz[i] <- mean(margyact[margypts==umyv[i]]) }
				margydf <- data.frame(x=margyposx,y=py,z=pz,w=pw)
			}
		}
	}
	if (logscale) {
		maindf$x <- exp(maindf$x)
		maindf$y <- exp(maindf$y)
		if (margins && exists("margxdf")) {
			margxdf$x <- exp(margxdf$x)
			margxdf$y <- exp(margxdf$y)
		}
		if (margins && exists("margydf")) {
			margydf$x <- exp(margydf$x)
			margydf$y <- exp(margydf$y)
		}
	}
	synp <- ggplot(maindf,aes_string(x="x",y="y",fill="z"))
	if (interpolate) {
		if (raster) { synp <- synp+geom_raster() }
		else { synp <- synp+geom_tile() }
		if (logscale) {
			if (margins) {
				if (exists("margxdf")) { synp <- synp+geom_tile(data=margxdf,height=margbarwidth) }
				if (exists("margydf")) { synp <- synp+geom_tile(data=margydf,width=margbarwidth) }
			}
			synp <- synp+scale_x_log10()+scale_y_log10()
		}
	} else {
		if (irreg) { synp <- synp+geom_polygon(aes_string(group="poly")) }
		else { synp <- synp+geom_tile() }
		if (logscale) {
			if (margins) {
				if (exists("margxdf")) { synp <- synp+geom_tile(aes_string(width="w"),data=margxdf,height=margbarwidth) }
				if (exists("margydf")) { synp <- synp+geom_tile(aes_string(height="w"),data=margydf,width=margbarwidth) }
			}
			synp <- synp+scale_x_log10()+scale_y_log10()
		}
	}
	return(synp)
}

formatResponseMap <- function(rmap,palette=NULL,cscale=NULL,xl=expression(Conc[A]),yl=expression(Conc[B]),zl="Activity") {
	fp <- rmap+labs(x=xl,y=yl)
	if (is.null(palette)) { palette <- c("#007F00", "green", "#7FFF00", "yellow", "#FF7F00", "red", "#7F0000") }
	if (!is.null(cscale)) {
		if (length(cscale)>1) { cscale <- cscale[1] }
		totmin <- min(rmap$data$z)
		totmax <- max(rmap$data$z)
		for (i in 1:length(rmap$layers)) {
			if (!is.null(rmap$layers[[i]]$data$z)) {
				totmin <- min(totmin,min(rmap$layers[[i]]$data$z))
				totmax <- max(totmax,max(rmap$layers[[i]]$data$z))
			}
		}
		dev <- max(totmax-cscale,cscale-totmin)+0.01
		clims <- c(cscale-dev,cscale+dev)
		fp <- fp + scale_fill_gradientn(zl,limits=clims,colours=palette)
	} else { fp <- fp + scale_fill_gradientn(zl,colours=palette) }
	
	return(fp)
}
formatDifferenceMap <- function(rmap,zcenter=NULL,xl=expression(Conc[A]),yl=expression(Conc[B]),zl="Diff") {
	fp <- rmap+labs(x=xl,y=yl)
	if (!is.null(zcenter)) { fp <- fp+scale_fill_gradient2(zl,midpoint=zcenter) }
	else { fp <- fp+scale_fill_gradient2(zl) }
	return(fp)
}

subdivide <- function(vect,ltarg=100) {
	k <- round(log2(ltarg/(length(unique(vect))-1)))
	dx <- (length(unique(vect))-1)*2^k
	xs <- seq(min(vect),max(vect),length=2*dx+1)
	xs <- xs[seq(2,length(xs)-1,by=2)]
	return(xs)
}
gaussInterp2d <- function(x,y,z,nx,ny,sigma=NULL) {
	if (length(x)!=length(y)) { stop("Vectors of x- and y-coordinates muat have the same length.") }
	if (is.null(sigma)) {
		sigx <- (max(x)-min(x))/(2*(length(unique(x))-1))
		sigy <- (max(y)-min(y))/(2*(length(unique(y))-1))
	} else { sigx <- sigy <- sigma }
	nz <- rep(0,length(nx))
	for (i in 1:length(nz)) {
		gv <- exp(-((nx[i]-x)^2/(2*sigx^2)+(ny[i]-y)^2/(2*sigy^2)))
		nz[i] <- sum(gv*z)/sum(gv)
	}
	return(nz)
}
gaussInterp1d <- function(x,y,nx,sigma=NULL) {
	if (is.null(sigma)) { sig <- (max(x)-min(x))/(2*(length(unique(x))-1)) }
	else { sig <- sigma }
	ny <- rep(0,length(nx))
	for (i in 1:length(ny)) {
		gv <- exp(-(nx[i]-x)^2/(2*sig^2))
		ny[i] <- sum(gv*y)/sum(gv)
	}
	return(ny)
}
linearInterpIrreg1d <- function(x,y,nx) {
	uxv <- sort(unique(x))
	ny <- rep(NA,length(nx))
	nwt <- rep(0,length(nx))
	for (i in 1:(length(uxv)-1)) {
		rel <- nx>=uxv[i] & nx<=uxv[i+1]
		pv <- (nx[rel]-uxv[i])/(uxv[i+1]-uxv[i])
		pv <- cbind(1-pv,pv)
		ny[rel] <- 0
		nwt[rel] <- 0
		for (j in 1:2) {
			for (k in which(x==uxv[i+j-1])) {
				ny[rel] <- ny[rel]+y[k]*pv[,j]
				nwt[rel] <- nwt[rel]+pv[,j]
			}
		}
	}
	nwt[nwt==0] <- 1
	ny <- ny/nwt
	ny[nx>max(x)] <- mean(y[x==max(x)])
	ny[nx<min(x)] <- mean(y[x==min(x)])
	return(ny)
}
linearInterpIrreg2d <- function(x,y,z,nx,ny) {
	nz <- rep(NA,length(nx))
	nwt <- rep(0,length(nx))
	upts <- unique(cbind(x,y))
	delts <- delaunay(upts)
	for (i in 1:nrow(delts)) {
		cupts <- upts[delts[i,],]
		s1 <- sign((ny-cupts[1,2])*(cupts[2,1]-cupts[1,1])-(nx-cupts[1,1])*(cupts[2,2]-cupts[1,2]))
		s2 <- sign((ny-cupts[2,2])*(cupts[3,1]-cupts[2,1])-(nx-cupts[2,1])*(cupts[3,2]-cupts[2,2]))
		s3 <- sign((ny-cupts[3,2])*(cupts[1,1]-cupts[3,1])-(nx-cupts[3,1])*(cupts[1,2]-cupts[3,2]))
		rel <- (s1>=0)&(s2>=0)&(s3>=0)
		cuvecs <- cupts[c(2,3,1),]-cupts
		pden <- (cuvecs[2,1]*cuvecs[3,2]-cuvecs[3,1]*cuvecs[2,2])
		pv <- array(0,dim=c(length(which(rel)),3))
		pv[,1] <- 1-(cuvecs[2,1]*(cupts[1,2]-ny[rel])-(cupts[1,1]-nx[rel])*cuvecs[2,2])/pden
		pv[,2] <- 1-(cuvecs[3,1]*(cupts[2,2]-ny[rel])-(cupts[2,1]-nx[rel])*cuvecs[3,2])/pden
		pv[,3] <- 1-(cuvecs[1,1]*(cupts[3,2]-ny[rel])-(cupts[3,1]-nx[rel])*cuvecs[1,2])/pden
		nz[rel] <- 0
		nwt[rel] <-0
		for (j in 1:3) {
			for (k in which(x==upts[delts[i,j],1] & y==upts[delts[i,j],2])) {
				nz[rel] <- nz[rel]+z[k]*pv[,j]
				nwt[rel] <- nwt[rel]+pv[,j]
			}
		}
	}
	hull <- convexHull(delts)
	ctx <- (max(x)+min(x))/2
	cty <- (max(y)+min(y))/2
	hth <- atan2(upts[hull,2]-cty,upts[hull,1]-ctx)
	nth <- atan2(ny-cty,nx-ctx)
	hth <- c(hth,hth[1])
	hull <- c(hull,hull[1])
	for (i in 1:(length(hull)-1)) {
		if (hth[i]<hth[i+1]) { rel <- is.na(nz) & nth>=hth[i] & nth<hth[i+1] }
		else { rel <- is.na(nz) & (nth>=hth[i] | nth<hth[i+1] ) }
		cupts <- rbind(c(ctx,cty),upts[c(hull[i],hull[i+1]),])
		cuvecs <- cupts[c(2,3,1),]-cupts
		pv <- (cuvecs[3,1]*(ny[rel]-cty)-(nx[rel]-ctx)*cuvecs[3,2])
		pv <- pv/((nx[rel]-ctx)*cuvecs[2,2]-cuvecs[2,1]*(ny[rel]-cty))
		pv <- cbind(pv,1-pv)
		nz[rel] <- 0
		nwt[rel] <-0
		for (j in 1:2) {
			for (k in which(x==upts[hull[i+j-1],1] & y==upts[hull[i+j-1],2])) {
				nz[rel] <- nz[rel]+z[k]*pv[,j]
				nwt[rel] <- nwt[rel]+pv[,j]
			}
		}
	}
	nwt[nwt==0] <- 1
	nz <- nz/nwt
	return(nz)
}
