# The EMbC Package for R
#
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic Bartumeus <fbartu@ceab.csic.es>
#
# EMbC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
#
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.


# binClst auxiliary functions
# ---------------------------

getColors <- function(k,fam='RdYlBu'){
	gry <- brewer.pal(4,'Greys')
	if (k==1){
		return(gry[4])
		}
	else if (k<=4){
		return(c(brewer.pal(4,fam)[c(2,1,3,4)],gry[4]))
		}
	else if(k<=8){
		return(c(brewer.pal(8,fam),gry[4]))
		}
	else return(c(rainbow(k),gry[4]))
	}

parSet <- function(mtx=NULL,widths=NULL,heights=NULL,bg=NULL,oma=c(1,1,1,1),mar=c(3,3,0.5,0.5),mgp=c(1.5,0.4,0),cex.main=1.0,cex.lab=1.0,cex.axis=0.8){
	parDef <- par(no.readonly=TRUE)
	par(oma=oma)
	par(mar=mar)
	par(mgp=mgp)
	par(cex.main=cex.main)
	par(cex.lab=cex.lab)
	par(cex.axis=cex.axis)
	if (!is.null(mtx)){
		if (is.null(widths)) widths <- rep(1,ncol(mtx))
		if (is.null(heights)) heights <- rep(1,nrow(mtx))
		layout(mtx,widths=widths,heights=heights)
		}
	if (!is.null(bg)) par(bg=bg)
	return(parDef)
	}

# binClstPath auxiliary functions
# -------------------------------

earthR <- function() return(6378137)

spanTime <- function(pth){
	n <- nrow(pth)
	return(as.numeric(difftime(pth$dTm[2:n],pth$dTm[1:(n-1)],units="secs")))}

loxDst <- function(pth){
	n <- nrow(pth)
	dltLon <- (pth$lon[2:n]-pth$lon[1:(n-1)])*pi/180
	dltLat <- (pth$lat[2:n]-pth$lat[1:(n-1)])*pi/180
	dltRad <- numeric((n-1))
	for (i in 1:(n-1)) {
		if (dltLat[i] == 0){
			q <- cos(pth$lat[i]*pi/180)}
		else{
			dltPhi <- log(tan(pi/4+pth$lat[i+1]*pi/360))-log(tan(pi/4+pth$lat[i]*pi/360))
			q <- dltLat[i]/dltPhi
			}
		dltRad[i] = sqrt(dltLat[i]^2+(q*dltLon[i])^2)
		}
	return(dltRad*earthR())}

loxTht <-function(pth){
	n <- nrow(pth)
	dltLon <- (pth$lon[2:n]-pth$lon[1:(n-1)])
	dltLat <- (pth$lat[2:n]-pth$lat[1:(n-1)])
	lclTht <- numeric((n-1))
	for (i in 1:(n-1)){
		tngTht <- ifelse (dltLon[i]==0,0,
			log(tan(pi/4+pth$lat[i+1]*pi/360)/tan(pi/4+pth$lat[i]*pi/360))/(dltLon[i]*pi/180))
		if (dltLon[i]==0 && dltLat[i]==0) lclTht[i] <- 2*pi
		else if (dltLon[i] == 0) lclTht[i] <- ifelse (dltLat[i]>0,0,pi)
		else if (dltLat[i] == 0) lclTht[i] <- ifelse (dltLon[i]>0,pi/2,3*pi/2)
		else {
			lclTht[i] <- ifelse (dltLon[i]>0,atan2(tngTht,1),atan2(-tngTht,-1))
			if (dltLat[i] > 0) {
				lclTht[i] <- 5*pi/2-lclTht[i]
				while (lclTht[i] > 2*pi)lclTht[i] <- lclTht[i]-2*pi}
			else {
				lclTht[i] <- pi/2-lclTht[i]}
			}
		}
	return(lclTht)}

getSpeed <- function(bCP){
	return(apply(cbind(bCP@dst,bCP@spn),1,function(x) if (x[2]>0) x[1]/x[2] else 0))}

getTurns <- function(bCP){
	Z <- numeric((nrow(bCP@pth)-1))
	Z[1] <- 0
	for (i in 2:length(Z)) {
		if (bCP@hdg[i] == 2*pi || bCP@hdg[i-1] == 2*pi) Z[i] <- 0
		else {
			Z[i] <- abs(bCP@hdg[i]-bCP@hdg[i-1])
			if (Z[i] > pi) Z[i] <- 2*pi - Z[i]
			}
		}
	return(c(Z,0))}

stCertainty <- function(bCP){
	uRef <- median(bCP@spn)
	if (uRef<=60){
		Ttbl <- sort(table(bCP@spn))
		iRef <- length(Ttbl)
		while (iRef>1 && uRef<=60){
			uRef <- as.numeric(names(Ttbl[iRef]))
			iRef <- iRef -1
			}
		}
	U <- sapply(bCP@spn,function(x) if (x>0) min(uRef/x,1) else 0)
	return(cbind(U,U))}

bCPStd <- function(pth){
	names(pth)[1:3] <- c('dTm','lon','lat')
	return(pth)}

setMarkerSizes <- function(bCP,nMarkerSizeClasses,minMarkerRadius,maxMarkerRadius,logDurations=TRUE){
	if(logDurations)
		durationVariable = log(bCP@midPoints@data$duration)
	else
		durationVariable = bCP@midPoints@data$duration
	range = maxMarkerRadius - minMarkerRadius
	if(range<0) return("Error: minmarkerRadius must be smaller than maxMarkerRadius. Please use new values.")
	if(length(nMarkerSizeClasses==1)){
		if(nMarkerSizeClasses==0){
			durationScaleFactor = range/(max(durationVariable)-min(durationVariable))
			circleRadiiMeters = (minMarkerRadius+durationScaleFactor*(durationVariable- min(durationVariable)))
		} else if(nMarkerSizeClasses==1){
			circleRadiiMeters = rep(minMarkerRadius + (maxMarkerRadius-minMarkerRadius)/2, length(durationVariable))
		}
		else if(nMarkerSizeClasses>1){
			markerSizeClasses = minMarkerRadius + (0:(nMarkerSizeClasses-1))*range/(nMarkerSizeClasses-1)
			circleRadiiMeters = markerSizeClasses[unclass(cut(durationVariable, nMarkerSizeClasses))]
		} else if(nMarkerSizeClasses < 0) return("Error: nMarkerSizeClasses must be >= 0. Please specify a new value.")
	} else if(length(nMarkerSizeClasses>1)){
		if(min(nMarkerSizeClasses)>min(durationVariable))
			nMarkerSizeClasses = c(min(durationVariable), nMarkerSizeClasses)
		else if(max(nMarkerSizeClasses)<max(durationVariable))
			nMarkerSizeClasses = c(nMarkerSizeClasses,max(durationVariable))
		markerSizeClasses = minMarkerRadius + (0:(length(nMarkerSizeClasses)-2))*range/(length(nMarkerSizeClasses)-2)
		circleRadiiMeters = markerSizeClasses[unclass(cut(durationVariable, nMarkerSizeClasses))]
	}
	return(circleRadiiMeters)}

getSolarPos <- function(pth,scv){
	solP <- solarpos(cbind(pth$lon,pth$lat),pth$dTm,proj4string=CRS("+proj=longlat +datum=WGS84"))
	if (scv=='both') return(solP)
	else if (scv=='azimuth') return(solP[,1])
	else if (scv=='height') return(solP[,2])
	else if (scv=='rheight') return(sin(solP[,2]*pi/180))
	else if (scv=='rheight2') return(sin(sin(solP[,2]*pi/180)*pi/2))
	else if (scv=='rheight3') return(sin(sin(sin(solP[,2]*pi/180)*pi/2)*pi/2))
	}

# Auxiliar Non-export functions
# -----------------------------

# binClstpath auxiliary functions

formatSecs <- function(secs){
	hr <- floor(secs/(60*60))
	min <-  floor((secs - (hr*60*60))/60)
	sec <- round(secs - ((hr*60*60)+(min*60)),digits=2)
	return(paste(hr,min,sec,sep=':'))}

formatMeters <- function(meters){
	if(meters < 1000) return(paste(round(meters, 0), "m"))
	else return(paste(round(meters/1000, 2), "km"))}

# format Tht parameters to given length and decimals.
frmTht <- function(Tht,d,w){
	paste(lapply(1:length(Tht$M),function(m){
		paste(formatC(Tht$M[m],format='f',digits=d,width=w),formatC(sqrt(Tht$S[m,m]),format='f',digits=d,width=w),sep=" ")}),collapse=" ")
	}

# clustering binary labels.
getkLbls <- function(bC,kNmbrs=FALSE)
	return(lapply(1:bC@k, function(k){
		bk <- paste(as.integer(rev(intToBits(k-1)[1:bC@m])),collapse="")
		if (kNmbrs)	{
			if (bC@k<9) paste(k,'.',gsub('1','H',gsub('0','L',bk)),sep="")
			else paste(formatC(k,width=2,flag='0'),'.',gsub('1','H',gsub('0','L',bk)),sep="")
			}
		else gsub('1','H',gsub('0','L',bk))
		}))

# max/min scaler
maxminScale <- function(x)
	return((x-min(x))/(max(x)-min(x)))

# get subPth with selected clusters
getSubPth <- function(bCP,showClst){
	subPth <- which(bCP@A %in% showClst)
	bCP@pth <- bCP@pth[subPth,]
	bCP@X <- bCP@X[subPth,]
	bCP@spn <- bCP@spn[subPth]
	bCP@dst <- bCP@dst[subPth]
	bCP@A <- bCP@A[subPth]
	return(bCP)}

# compute proportional limits for path view
getPropLims <- function(pth, a, b){
	plims <- list(x=c(min(pth$lon[a:b]), max(pth$lon[a:b])),
	y=c(min(pth$lat[a:b]), max(pth$lat[a:b])))
	shorter <- which.min(abs(sapply(plims, diff)))
	larger <- which.max(abs(sapply(plims, diff)))
	while (abs(diff(plims[[shorter]])) < abs(diff(plims[[larger]]))) plims[[shorter]] <- plims[[shorter]] + c(-0.0001,+0.0001)
	return(plims)
}

# bivariate binary clustering scatterplot with reference lines and legend
sctr2D <- function(bC){
	if (length(bC@A)==0){
		plot(bC@X[,1],bC@X[,2],pch=20,xlab=colnames(bC@X)[1],ylab=colnames(bC@X)[2])}
	else {
		plot(bC@X[,1],bC@X[,2],col=bC@C[bC@A],pch=20,xlab=colnames(bC@X)[1],ylab=colnames(bC@X)[2])}
	lines(c(bC@R[1,3],bC@R[1,3]),c(bC@R[1,2],bC@R[1,4]),col='grey')
	lines(c(bC@R[1,1],bC@R[1,3]),c(bC@R[1,4],bC@R[1,4]),col='grey')
	lines(c(bC@R[4,1],bC@R[4,1]),c(bC@R[4,2],bC@R[4,4]),col='grey')
	lines(c(bC@R[4,1],bC@R[4,3]),c(bC@R[4,2],bC@R[4,2]),col='grey')
	legend("topright",legend=c(getkLbls(bC),'NC'),col=bC@C,cex=0.8,lwd=3,text.font=1,bty='n')
	}


# multivariate binary clustering scatterplot
sctr3D <- function(obj,showVars=numeric(),showClst=numeric()){
	# set vars and clusters subset
	lims <- apply(obj@X,2,range)
	if (length(showVars)>=3) m <- showVars[1:3]
	else m <- c(1,2,3)
	X <- obj@X
	A <- obj@A
	if (length(showClst)>0) {
		X <- X[which(A %in% showClst),]
		A <- A[which(A %in% showClst)]}
	# titles and axes labels
	par(mar=c(4,4,2,0.5))
	mttl <- paste(colnames(obj@X)[m[1]],c(":LOW",":HIGH"),sep='')
	labs <- paste(colnames(obj@X)[m[2:3]],sep='')
	# plot scatterplpot for LOW  values of the splitting variable
	klow <- which(substr(getkLbls(obj),m[1],m[1])=='L')
	Xlow <- X[which(A %in% klow),]
	plot(Xlow[,m[2]],Xlow[,m[3]],col=obj@C[A[which(A %in% klow)]],pch=20,xlim=lims[,m[2]],ylim=lims[,m[3]],main=mttl[1],xlab=labs[1],ylab=labs[2],cex.main=1.0,font.main=1)
	# plot scatterplpot for HIGH values of the splitting variable
	khgh <- which(substr(getkLbls(obj),m[1],m[1])=='H')
	Xhgh <- X[which(A %in% khgh),]
	plot(Xhgh[,m[2]],Xhgh[,m[3]],col=obj@C[A[which(A %in% khgh)]],pch=20,xlim=lims[,m[2]],ylim=lims[,m[3]],main=mttl[2],xlab=labs[1],ylab=labs[2],cex.main=1.0,font.main=1)
	# plot legend
	par(mar=c(4,1,0.5,0.5))
	plot(0,axes=FALSE,xlab="",ylab="",pch='')
	obj@C[which(lapply(1:obj@k,function(k) length(obj@A[which(obj@A==k)])==0)==TRUE)] <- brewer.pal(4,'Greys')[3]
	if (length(showClst)>0) obj@C[which(lapply(1:obj@k,function(k) !(k%in%showClst))==TRUE)] <- brewer.pal(4,'Greys')[3]
	legend("center",legend=c(getkLbls(obj,kNmbrs=TRUE),'NC'),col=obj@C,cex=0.8,lwd=3,text.font=1,bty='n')
	}
