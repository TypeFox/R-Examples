
dendro = function(x,k = 30, color.id = k-2, label = FALSE, opts = list(), min.gap = 0.01,spline = FALSE, ...){
	stopifnot( !is.null(x$height) )
	stopifnot( length(x$height) >= k )
	stopifnot(color.id <= k)
	
	if(!"alpha" %in% names(opts)){
		opts$alpha = 1	
	}
	if(!"lwd" %in% names(opts)){
		opts$lwd= 1	
	}
	if(!"ps" %in% names(opts)){
		opts$ps = 1	
	}

	
	xcoords <- rev(tail(x$height,k-1))
	xcoords <- xcoords/max(x$height) 
	
	if(min.gap > 0){
		# combine all joins with gaps <= min.gap
		dx <- diff(xcoords)
		if(any(dx < min.gap)){
			dx <- dist(xcoords)
			hdx <- hclust(dx,method="single")
			#ids <- subtree(hdx, h = min.gap)$data
			ids <- cutree(hdx, h = min.gap)
			k.vec <- unlist(tapply(2:k,ids, function(s){
				if(length(s) > 1){
					return(c(s[1], s[1]))	
				}
				return(s)
			}))
#			ycoords <-	by(t(ycoords),ids,function(s){
#				ds <- dim(s)
#				s2 <- s[1,]
#				if(ds[1]>1){
#					s2 <- rbind(s2,s[1,])
#				}
#					return(s2)	
#			},simplify=FALSE)
#			ycoords <- t(do.call(rbind,ycoords))
		}
		xcoords <- sort(unlist(tapply(xcoords,ids, function(v) unique(range(v)))),decreasing=TRUE)
		k <- length(k.vec)+1
		color.id <- min(k-1,color.id)
	}else{
		ids <- 1:(k-1)
		k.vec <- 2:k	
	}
#	cs <- NULL
#for(i in 2:k) cs <- cbind(cs, subtree(x,k=i)$data)
cs <- sapply(k.vec, function(i){#
	as.integer(subtree(x,k=i)$data)
})
	cs <- as.data.frame(cs)
	names(cs) <- paste("V",2:k,sep="")
	ss <- scpcp(cs,sort.individual=TRUE, plot = FALSE)
	
	
	
	ycoords <- data.matrix(ss[,1:k])
	ycoords <- apply(ycoords,2,function(s){
		y <- s - min(s)
		y <- y/max(y)
		return(y)
	})
	
	
	
	
	# components with either 1 or two y-coordinates.
	# two occurs if 2 or more coordinates have been collapsed 
	comp <- 1+(count(ids)[,2]>1)
	
	ycoords <- cbind(ycoords, as.integer(as.matrix(ss[,color.id])))
	colv <- alpha(rainbow_hcl(color.id+1),opts$alpha)
	#dev.new()
	grid.newpage()
	vp <- viewport( 0.5, 0.5, width = 0.9, height = 0.9)
	pushViewport(vp)

	apply(ycoords,1,function(s){
		grid.points(x=s[-k], y=xcoords, pch=19, gp=gpar(col = colv[s[k]]),size=unit(0.002*opts$ps,"npc"))
		
		if(spline){
			psp <- pairspline(xcoords,s[-k])
			grid.lines(x=psp[,2], y=psp[,1], gp=gpar(col = colv[s[k]],lwd=opts$lwd))
		}else{
			grid.lines(x=s[-k], y=xcoords, gp=gpar(col = colv[s[k]],lwd=opts$lwd))
		}
		
		
	})
	
	j = 0
	for(i in seq_along(comp)){
		j <- j+comp[i]
		grp <- ss[,3*k-3+j]
		ht <- xcoords[j]
	
		if(comp[i] > 1){
			crd <- cbind(ycoords[,j],ycoords[,j-1])
			ht2 <- xcoords[j-1]
			by(crd,grp,function(w){
				grid.rect( x = min(w[,1]), y = ht, width= max(w[,2])-min(w[,1]),height=ht2-ht,
				 just=c("left","bottom"),gp=gpar(fill=alpha(grey(0.2),0.1)))#gp=gpar(fill=alpha(grey(0.2),0.05))	
			})
		}else{
			crd <- ycoords[,j]
			tapply(crd,grp,function(w){
				grid.lines( x = c(min(w),max(w)), y = c(ht,ht),gp = gpar(lwd=2))	
			})
		}
	}
		
	popViewport()
	
	vpy <- viewport( 0.025, 0.5, width = 0.05, height = 0.9)
	pushViewport(vpy)
	unreg.axis(x0=0.9, y0=0, ticks = round(c(0,xcoords),2), rot=90)
	popViewport()
	
	if( label ){
			pushViewport(vp)
			if(is.null(x$order)){
				x$order <- 1:length(x$height)	
			}
			if(is.null(x$labels)){
				x$labels <- 1:length(x$height)	
			}
			
			nmz <- x$labels[x$order]
			grid.text(label = nmz, x=ycoords[,k-1], y=rep(0.9*min(xcoords),length(x$height)),just="right",rot=90, gp=gpar(col=colv[ycoords[,k]]))
			popViewport()
	}
}



unreg.axis = function(x0,y0,rot, ticks, len = 1, ltm = 1/30, clockwise = FALSE, lab.cex = 1.0, keep=7, col.axis = 1){
		
	grid.lines(x=c(x0, x0+len*cos(rot/180*pi)), y = c(y0, y0+len*sin(rot/180*pi)),gp=gpar(col=col.axis))
	n <- length(ticks)-1
	
	sgn <- ifelse(clockwise,-1,1)
	
	xs <- ticks - min(ticks)
	xs <- xs/max(xs)
	
	x.ticks.1 <- x0+xs*len*cos(rot/180*pi)
	y.ticks.1 <- y0+xs*len*sin(rot/180*pi)
	
	x.ticks.2 <- x.ticks.1 + sgn*ltm*len*cos(rot/180*pi+pi/2)
	y.ticks.2 <- y.ticks.1 + sgn*ltm*len*sin(rot/180*pi+pi/2)
	
	mapply( function(x0,y0,x1,y1){
		grid.lines(c(x0,x1),c(y0,y1),gp=gpar(col=col.axis))
		},x0 = x.ticks.1, x1 = x.ticks.2, y0= y.ticks.1, y1 = y.ticks.2)
		if(rot <= 135 & rot > 45) just = ifelse(clockwise,"left","right")
		if(rot <= 225 & rot > 135) just = ifelse(clockwise,"bottom","top")
		if(rot <= 315 & rot > 225) just = ifelse(clockwise,"right","left")
		if(rot <= 45 | rot > 315) just = ifelse(clockwise,"left","right")#top/bottom
	
#	ii <- round( (n-1)/4 )
#	if(ii > 0){
#		ids <- which( floor(1:(n-2)+ii/2) %% ii != 0)
#		ticks[ids+1] = ""	
#	}
	
#	tmp <- n / c(5,6)
#	tick.step <- trunc(tmp)[which.min( tmp - trunc(tmp)  )]
if(n > keep){
	keep <- round(seq(1,n,n/(keep-1)))
	ticks[-c(1,n+1,keep)] = ""
}		
	grid.text(ticks, x.ticks.2, y.ticks.2, just = just, rot = 0, gp=gpar(fontsize=lab.cex*10,col=col.axis))
}



pairspline = function(x,y,n=20,e = 12){
	
	n <- length(x)
	cc<-sapply(2:n,function(w){
		x2 <- x[w]
		x1 <- x[w-1]
		y2 <- y[w]
		y1 <- y[w-1]
		dx <- x2-x1
		dy <- y2-y1
		xc <- x1 + seq(0,1,1/n)*dx
		yc <- y1 + (1/(1+exp(-e*seq(-0.5,0.5,1/n))))*dy
		return(rbind(xc,yc))
	},simplify=FALSE)
cc <- do.call(cbind,cc)


return(t(cc))
}

