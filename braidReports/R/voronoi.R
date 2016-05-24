voronoi <- function(pts,bbox=NULL) {
	if (is.null(bbox)) {
		rx <- range(pts[,1])
		ry <- range(pts[,2])
		f <- 1/(2*sqrt(nrow(pts)))
		bbox <- c((1+f)*rx-f*rx[c(2,1)],(1+f)*ry-f*ry[c(2,1)])
	}
	
	tris <- delaunay(pts,bbox)
	delsegs <- array(0,dim=c(0,7))
	for (ti in 1:nrow(tris)) {
		v1 <- pts[tris[ti,2],]-pts[tris[ti,1],]
		v2 <- pts[tris[ti,3],]-pts[tris[ti,1],]
		D <- 2*(v1[1]*v2[2]-v1[2]*v2[1])
		if (D==0) { next }
		cx <- pts[tris[ti,1],1]+(v2[2]*sum(v1^2)-v1[2]*sum(v2^2))/D
		cy <- pts[tris[ti,1],2]+(v1[1]*sum(v2^2)-v2[1]*sum(v1^2))/D
		for (e in 1:3) {
			eseg <- c(tris[ti,c(e,(e%%3)+1)],Inf,Inf,cx,cy,0)
			if (nrow(delsegs)>0) { match <- which(delsegs[,1]==eseg[2]&delsegs[,2]==eseg[1]) }
			if (nrow(delsegs)==0 || length(match)==0) {
				delsegs <- rbind(delsegs,eseg)
			} else if (nrow(delsegs)>0) { delsegs[match,3:4] <- c(cx,cy) }
		}
	}
	
	adelsegs <- delsegs
	outs <- is.infinite(adelsegs[,3])
	adelsegs[,7] <- 1*outs
	adelsegs[outs,3:4] <- (pts[adelsegs[outs,1],]+pts[adelsegs[outs,2],])/2
	outsign <- sign((adelsegs[outs,6]-adelsegs[outs,4])*(pts[adelsegs[outs,2],1]-pts[adelsegs[outs,1],1])-
						(adelsegs[outs,5]-adelsegs[outs,3])*(pts[adelsegs[outs,2],2]-pts[adelsegs[outs,1],2]))
	coinc <- which(outs)[outsign==0]
	adelsegs[coinc,3:4] <- adelsegs[coinc,3:4]+cbind(pts[adelsegs[coinc,2],2]-pts[adelsegs[coinc,1],2],
													pts[adelsegs[coinc,1],1]-pts[adelsegs[coinc,2],1])
	adelsegs[which(outs)[outsign<0],7] <- -1

	inind <- (adelsegs[,7]!=0|adelsegs[,3]<bbox[1]|adelsegs[,3]>bbox[2]|adelsegs[,4]<bbox[3]|adelsegs[,4]>bbox[4])
	inpts <- adelsegs[inind,]
	ibx1 <- bbox[1*xor(inpts[,5]<=inpts[,3],inpts[,7]<0)+1]
	ibt1 <- (ibx1-inpts[,3])/(inpts[,5]-inpts[,3])
	ibt1[which(inpts[,5]==inpts[,3])] <- Inf
	iby1 <- inpts[,4]+ibt1*(inpts[,6]-inpts[,4])
	val1 <- (iby1>=bbox[3] & iby1<=bbox[4] & (inpts[,7]>0|ibt1>=0) & ((inpts[,7]<0&ibt1>1)|(inpts[,7]>=0&ibt1<=1)))
	inpts[val1,3:4] <- cbind(ibx1[val1],iby1[val1])
	adelsegs <- rbind(adelsegs[!inind,],inpts[val1,])
	inpts <- inpts[!val1,]
	if (length(inpts)>0) {
		if (length(inpts)==7) { inpts <- array(inpts,dim=c(1,7)) }
		iby2 <- bbox[1*xor(inpts[,6]<=inpts[,4],inpts[,7]<0)+3]
		ibt2 <- (iby2-inpts[,4])/(inpts[,6]-inpts[,4])
		ibt2[which(inpts[,6]==inpts[,4])] <- Inf
		ibx2 <- inpts[,3]+ibt2*(inpts[,5]-inpts[,3])
		val2 <- (ibx2>=bbox[1] & ibx2<=bbox[2] & (inpts[,7]>0|ibt2>=0) & ((inpts[,7]<0&ibt2>1)|(inpts[,7]>=0&ibt2<=1)))
		inpts[val2,3:4] <- cbind(ibx2[val2],iby2[val2])
		adelsegs <- rbind(adelsegs,inpts[val2,])
	}

	otind <- (adelsegs[,5]<bbox[1]|adelsegs[,5]>bbox[2]|adelsegs[,6]<bbox[3]|adelsegs[,6]>bbox[4])
	otpts <- adelsegs[otind,]
	if (length(otpts)>0) {
		if (length(otpts)==7) { otpts <- array(otpts,dim=c(1,7)) }
		ibx1 <- bbox[1*(otpts[,5]>=otpts[,3])+1]
		ibt1 <- (ibx1-otpts[,3])/(otpts[,5]-otpts[,3])
		ibt1[otpts[,5]==otpts[,3]] <- Inf
		iby1 <- otpts[,4]+ibt1*(otpts[,6]-otpts[,4])
		val1 <- (iby1>=bbox[3] & iby1<=bbox[4] & ibt1>=0 & ibt1<=1)
		otpts[val1,5:6] <- cbind(ibx1[val1],iby1[val1])
		adelsegs <- rbind(adelsegs[!otind,],otpts[val1,])
		otpts <- otpts[!val1,]
		if (length(otpts)>0) {
			if (length(otpts)==7) { otpts <- array(otpts,dim=c(1,7)) }
			iby2 <- bbox[1*(otpts[,6]>=otpts[,4])+3]
			ibx2 <- otpts[,3]+(iby2-otpts[,4])*(otpts[,5]-otpts[,3])/(otpts[,6]-otpts[,4])
			otpts[,5:6] <- cbind(ibx2,iby2)
			adelsegs <- rbind(adelsegs,otpts[,])
		}
	}
	adelsegs <- adelsegs[adelsegs[,3]!=adelsegs[,5]|adelsegs[,4]!=adelsegs[,6],]

	allpolys <- array(0,dim=c(0,3))
	for (p in 1:nrow(pts)) {
		psegs <- rbind(adelsegs[adelsegs[,1]==p,2:6],adelsegs[adelsegs[,2]==p,c(1,5,6,3,4)])
		if (length(psegs)==5) { psegs <- array(psegs[c(2:5)],dim=c(1,4)) }
		else {
			spks <- pts[psegs[,1],]
			spks[,1] <- spks[,1]-pts[p,1]
			spks[,2] <- spks[,2]-pts[p,2]
			psegs <- psegs[order(atan2(spks[,2],spks[,1])),2:5]
		}
		psegs <- rbind(psegs,psegs[1,])
		ppoly <- array(0,dim=c(0,2))
		for (s in 1:nrow(psegs)) {
			if (s>1 && (psegs[s,1]!=psegs[s-1,3] || psegs[s,2]!=psegs[s-1,4])) {
				ibx <- psegs[s-1,3]
				iby <- psegs[s-1,4]
				for (iter in 1:4) {
					if (ibx==psegs[s,1] || iby==psegs[s,2]) { break }
					if (ibx==bbox[1]&&iby!=bbox[3]) { iby <- bbox[3] }
					else if (iby==bbox[3]&&ibx!=bbox[2]) { ibx <- bbox[2] }
					else if (ibx==bbox[2]&&iby!=bbox[4]) { iby <- bbox[4] }
					else { ibx <- bbox[1] }
					ppoly <- rbind(ppoly,c(ibx,iby))
				}
				ppoly <- rbind(ppoly,psegs[s,1:2])
			}
			if (s<nrow(psegs)) { ppoly <- rbind(ppoly,psegs[s,3:4]) }
		}
		allpolys <- rbind(allpolys,cbind(ppoly,rep(p,times=nrow(ppoly))))
	}
	allpolys <- as.data.frame(allpolys)
	names(allpolys) <- c("x","y","poly")
	return(allpolys)
}

delaunay <- function(pts,bbox=NULL) {
	if (is.null(bbox)) {
		rx <- range(pts[,1])
		ry <- range(pts[,2])
		f <- 1/(2*sqrt(nrow(pts)))
		bbox <- c((1+f)*rx-f*rx[c(2,1)],(1+f)*ry-f*ry[c(2,1)])
	}

	rord <- 1:nrow(pts)
	# rord <- sample(nrow(pts))
	pts <- pts[rord,]
	
	p3 <- which((pts[3:nrow(pts),1]-pts[1,1])*(pts[2,2]-pts[1,2])-(pts[3:nrow(pts),2]-pts[1,2])*(pts[2,1]-pts[1,1])!=0)
	if (length(p3)==0) { stop("All points colinear.  This is not handled yet.") }
	else { p3 <- p3[1]+2 }
	if (p3!=3) {
		pts[c(3,p3),] <- pts[c(p3,3),]
		rord[c(3,p3)] <- rord[c(p3,3)]
	}
	
	tsign <- det(as.matrix(cbind(pts[1:3,],rep(1,times=3))))
	if (tsign>0) {
		tris <- array(c(1,2,3),dim=c(1,3))
		hull <- array(c(1,2,3,2,3,1),dim=c(3,2))
	} else {
		tris <- array(c(1,3,2),dim=c(1,3))
		hull <- array(c(1,3,2,3,2,1),dim=c(3,2))
	}

	# Add points iteratively
	# For each point
	eps <- (10^-6)*min(bbox[4]-bbox[3],bbox[2]-bbox[1])
	for (p in 4:nrow(pts)) {
		hullsign <- (pts[p,1]-pts[hull[,1],1])*(pts[hull[,2],2]-pts[hull[,1],2])-
						(pts[p,2]-pts[hull[,1],2])*(pts[hull[,2],1]-pts[hull[,1],1])
		if (length(which(hullsign>0))==0) {
			if (length(which(hullsign==0))>0) {
				# If point lies ON a hull edge, form two new triangles and hull edges
				if (length(which(hullsign==0))>1) {
					eproj <- (pts[p,1]-pts[hull[,1],1])*(pts[hull[,2],1]-pts[hull[,1],1])+
								(pts[p,2]-pts[hull[,1],2])*(pts[hull[,2],2]-pts[hull[,1],2])
					edot <-(pts[hull[,2],1]-pts[hull[,1],1])*(pts[hull[,2],1]-pts[hull[,1],1])+
								(pts[hull[,2],2]-pts[hull[,1],2])*(pts[hull[,2],2]-pts[hull[,1],2])
					wh <- which(hullsign==0 & eproj>0 & eproj<edot)[1]
				} else { wh <- which(hullsign==0)[1] }
				itri <- which((tris[,1]==hull[wh,1] & tris[,2]==hull[wh,2])|
								(tris[,2]==hull[wh,1] & tris[,3]==hull[wh,2])|
								(tris[,3]==hull[wh,1] & tris[,1]==hull[wh,2]))
				tsd <- which(tris[itri,]==hull[wh,1])
				fliplist <- rbind(c(tris[itri,c(tsd%%3+1,(tsd+1)%%3+1)],p),
									c(tris[itri,c((tsd+1)%%3+1,(tsd+2)%%3+1)],p))
				if (wh==1) { hull <- rbind(c(hull[1,1],p),c(p,hull[1,2]),hull[-1,]) }
				else if (wh==nrow(hull)) { hull <- rbind(hull[-nrow(hull),],c(hull[nrow(hull),1],p),c(p,hull[nrow(hull),2])) }
				else { hull <- rbind(hull[1:(wh-1),],c(hull[wh,1],p),c(p,hull[wh,2]),hull[(wh+1):nrow(hull),]) }
			} else {
				# If point is in convex hull/inside an existing triangle
				# Add three internal triangles to triangle list
				# Add three edges of existing triangle to "flip list"
				if (nrow(tris)==1) { tripts <- array(c(pts[tris[,1],],pts[tris[,2],],pts[tris[,3],]),dim=c(1,6)) }
				else { tripts <- cbind(pts[tris[,1],],pts[tris[,2],],pts[tris[,3],]) }
				pvec <- c(-pts[p,2],pts[p,1],pts[p,2],-pts[p,1],0,0)
				trisign <- tripts %*% cbind(pvec,pvec[c(5,6,1,2,3,4)],pvec[c(3,4,5,6,1,2)])
				trisign <- trisign+tripts[,c(1,3,5)]*tripts[,c(4,6,2)]
				trisign <- trisign-tripts[,c(3,5,1)]*tripts[,c(2,4,6)]
				itri <- which(trisign[,1]>=-eps & trisign[,2]>=-eps & trisign[,3]>=-eps)[1]
				fliplist <- rbind(c(tris[itri,1:2],p),c(tris[itri,2:3],p),c(tris[itri,c(3,1)],p))
			}
			tris <- rbind(tris[-itri,],fliplist)
			for (h in 1:nrow(hull)) {
				if (length(fliplist)==0) { break }
				for (f in 1:nrow(fliplist)) {
					if (fliplist[f,1]==hull[h,1] && fliplist[f,2]==hull[h,2]) {
						fliplist <- fliplist[-f,]
						if (length(fliplist)==3) { fliplist <- array(fliplist,dim=c(1,3)) }
						break
					}
				}
			}
		} else {
			# If point is outside convex hull,
			# Determine convex hull edges that face point
			# Add triangle with new point for each hull edge
			# Add previous hull edges to "flip list"
			
			if (hullsign[1]<=0 && hullsign[nrow(hull)]<=0) {
				fst <- min(which(hullsign>0))
				lst <- max(which(hullsign>0))
				sg1 <- c(hull[fst-1,2],p)
				sg2 <- c(p,hull[lst+1,1])
				nsegs <- hull[fst:lst,]
				hull <- rbind(hull[1:(fst-1),],sg1,sg2,hull[(lst+1):nrow(hull),])
			} else {
				fst <- min(which(hullsign<=0))
				lst <- max(which(hullsign<=0))
				sg1 <- c(hull[lst,2],p)
				if (lst<nrow(hull)) { nsegs <- hull[(lst+1):nrow(hull),] }
				else { nsegs <- array(0,dim=c(0,2)) }
				sg2 <- c(p,hull[fst,1])
				if (fst>1) { nsegs <- rbind(nsegs,hull[1:(fst-1),]) }
				hull <- rbind(hull[fst:lst,],sg1,sg2)
			}
			if (length(nsegs)==2) { fliplist <- c(nsegs[c(2,1)],p) }
			else if (length(nsegs)>2) { fliplist <- cbind(nsegs[,c(2,1)],rep(p,times=nrow(nsegs))) }
			tris <- rbind(tris,fliplist)
		}
		if (length(fliplist)==0) { next }
		else if (length(fliplist)==3) { fliplist <- array(fliplist,dim=c(1,3)) }
		
		
		# Loop through flip list
		iter <- 1
		while (iter<=nrow(fliplist)) {
			# Check if edge needs to be flipped
			# If edge needs to be flipped, remove existing triangle pair
			# Add new triangle pair
			# Add four outside edges of triangle pair to flip list
			fseg <- fliplist[iter,]
			tri2 <- which(tris[,2]==fseg[1]&tris[,1]==fseg[2])
			if (length(tri2)>0) { fseg[4] <- tris[tri2[1],3] }
			else {
				tri2 <- which(tris[,3]==fseg[1]&tris[,2]==fseg[2])
				if (length(tri2)>0) { fseg[4] <- tris[tri2[1],1] }
				else {
					tri2 <- which(tris[,1]==fseg[1]&tris[,3]==fseg[2])
					fseg[4] <- tris[tri2[1],2]
				}
			}
			tri2 <- tri2[1]
			fmat <- cbind(pts[fseg[1:3],1]-pts[fseg[4],1],pts[fseg[1:3],2]-pts[fseg[4],2],
							pts[fseg[1:3],1]^2-pts[fseg[4],1]^2+pts[fseg[1:3],2]^2-pts[fseg[4],2]^2)
			
			if (det(as.matrix(fmat))>0) {
				tri1 <- which(tris[,1]==fseg[1]&tris[,2]==fseg[2])
				if (length(tri1)==0) {
					tri1 <- which(tris[,2]==fseg[1]&tris[,3]==fseg[2])
					if (length(tri1)==0) { tri1 <- which(tris[,3]==fseg[1]&tris[,1]==fseg[2]) }
				}
				tri1 <- tri1[1]
				tris <- rbind(tris[-c(tri1,tri2),],fseg[c(2,3,4)],fseg[c(1,4,3)])
				nsegs <- array(fseg[c(2,3,1,4,3,1,4,2,4,4,3,3)],dim=c(4,3))
				for (h in 1:nrow(hull)) {
					if (length(nsegs)==0) { break }
					for (n in 1:nrow(nsegs)) {
						if (nsegs[n,1]==hull[h,1] && nsegs[n,2]==hull[h,2]) {
							nsegs <- nsegs[-n,]
							if (length(nsegs)==3) { nsegs <- array(nsegs,dim=c(1,3)) }
							break
						}
					}
				}
				if (iter<nrow(fliplist)) {
					fliplist <- fliplist[(iter+1):nrow(fliplist),]
					if (length(fliplist)==3) { fliplist <- array(fliplist,dim=c(1,3)) }
					for (f in 1:nrow(fliplist)) {
						if (length(nsegs)==0) { break }
						for (n in 1:nrow(nsegs)) {
							if ((nsegs[n,1]==fliplist[f,1] && nsegs[n,2]==fliplist[f,2]) ||
									(nsegs[n,1]==fliplist[f,2] && nsegs[n,2]==fliplist[f,1])) {
								nsegs <- nsegs[-n,]
								if (length(nsegs)==3) { nsegs <- array(nsegs,dim=c(1,3)) }
								break
							}
						}
					}
					fliplist <- rbind(fliplist,nsegs)
				} else { fliplist <- nsegs }
				iter <- 1
			} else { iter <- iter+1 }
		}
	}
	
	tris <- array(rord[tris],dim=dim(tris))
	return(tris)
}

convexHull <- function(tris) {
	hull <- c()
	for (i in 1:nrow(tris)) {
		for (j in 1:3) {
			ind <- which((tris[,1]==tris[i,(j%%3)+1] & tris[,2]==tris[i,j])|(tris[,2]==tris[i,(j%%3)+1]
						& tris[,3]==tris[i,j])|(tris[,3]==tris[i,(j%%3)+1] & tris[,1]==tris[i,j]))
			if (length(ind)==0) {
				hull <- c(tris[i,j],tris[i,(j%%3)+1])
				break
			}
		}
		if (length(hull)>0) { break }
	}
	for (i in 1:max(as.vector(tris))) {
		ch <- hull[length(hull)]
		posts <- c(tris[tris[,1]==ch,2],tris[tris[,2]==ch,3],tris[tris[,3]==ch,1])
		pres <- c(tris[tris[,1]==ch,3],tris[tris[,2]==ch,1],tris[tris[,3]==ch,2])
		nh <- setdiff(posts,pres)[1]
		if (nh==hull[1]) { break }
		else { hull <- c(hull,nh) }
	}
	return(hull)
}
