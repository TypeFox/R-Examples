plotDensity3d.fnc<-function(x,
			y,
			plot.type="contour",
			color="terrain",
			xlab=NULL, 
			ylab=NULL, 
			zlab=NULL, 
			main=NULL, 
			cex=1,
			alpha=1,
			lit=TRUE,
			theta=0,
			phi=0,
			bw="nrd0",
			adjust=1,
			kernel=c("gaussian","epanechnikov","rectangular","triangular","biweight","cosine","optcosine"),
			weights=NULL,
			window=kernel,
			width,
			give.Rkern=FALSE,
			n=50,
			from,
			to,
			cut=3,
			na.rm=FALSE,
                        ...
){
	# get unique values
	x<-sort(unique(x))
	y<-sort(unique(y))
	# get densities
	xd<-density(x=x,bw=bw,adjust=adjust,kernel=kernel,weights=weights,
		window=kernel,width=width,give.Rkern=give.Rkern,n=n,from=from,
		to=to,cut=cut,na.rm=na.rm)
	yd<-density(x=y,bw=bw,adjust=adjust,kernel=kernel,weights=weights,
		window=kernel,width=width,give.Rkern=give.Rkern,n=n,from=from,
		to=to,cut=cut,na.rm=na.rm)
	# get x*y matrix
	mat<-outer(xd$y,yd$y)

	# set labels if NULL
	if(is.null(xlab)){
		xlab=paste("x: N =",xd$n,"Bandwidth =",round(xd$bw,4),sep=" ")
	}

	if(is.null(ylab)){
		ylab=paste("y: N =",yd$n,"Bandwidth =",round(yd$bw,4),sep=" ")
	}

	if(is.null(zlab)){
		zlab<-"Density"
	}

	if(is.null(main)){
		if(plot.type=="contour"){
			main=zlab
		}else{
			main=""
		}
	}

	# contour plot
	if(plot.type=="contour"){
		#contourlevels = seq(zlim[1], zlim[2], by=contourstepsize)
			
		# Determine color.
        	if(color=="heat"){
            		pal=heat.colors(50)
            		con.col=3
        	}else if(color=="topo"){
            		pal=topo.colors(50)
            		con.col=2
        	}else if(color=="cm"){
            		pal=cm.colors(50)
            		con.col=1
        	}else if(color=="terrain"){
            		pal=terrain.colors(50)
            		con.col=2
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		pal=gray(seq(0.1,0.9,length=50))
            		con.col=1
        	}else{
			stop("please specify one of the following colors: heat, topo, cm, terrain, gray, or bw")
		}
		
		# plot
		image(x=seq(min(x),max(x),length=nrow(mat)),y=seq(min(y),max(y),length=ncol(mat)),
			z=mat,col=pal,main=main,xlab=xlab,ylab=ylab,axes=TRUE,cex.main=cex,
			cex.lab=cex,cex.axis=cex,...)
		contour(x=seq(min(x),max(x),length=nrow(mat)),y=seq(min(y),max(y),length=ncol(mat)),
			z=mat,add=TRUE,axes=FALSE,...)
		box()

		# return info
                return(invisible(list(x=x,y=y,xd=xd,yd=yd,mat=mat,col=pal)))
	}

	# perspective plot
	if(plot.type=="persp"){
		# the color portion of this code is adapted from the persp() help page
		#par(bg="white")
		nrz<-nrow(mat)
		ncz<-ncol(mat)

		# Create a function interpolating colors in the range of specified colors
        	if(color=="heat"){
            		jet.colors<-colorRampPalette(heat.colors(50))
        	}else if(color=="topo"){
			#jet.colors <- colorRampPalette( c("purple","blue", "green","yellow","red","white") ) 
			jet.colors <- colorRampPalette(topo.colors(50)) 
        	}else if(color=="cm"){
            		jet.colors<-colorRampPalette(cm.colors(50))
        	}else if(color=="terrain"){
            		jet.colors<-colorRampPalette(terrain.colors(50))
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		jet.colors<-colorRampPalette(gray(seq(0.1,0.9,length=7)))
        	}

		# Generate the desired number of colors from this palette
		nbcol<-100
		color<-jet.colors(nbcol)

		# Compute the z-value at the facet centres
		zfacet<-mat[-1,-1]+mat[-1,-ncz]+mat[-nrz,-1]+mat[-nrz,-ncz]

		# Recode facet z-values into color indices
		facetcol<-cut(zfacet,nbcol)

		# plot
		persp(x=seq(min(x),max(x),length=nrow(mat)),y=seq(min(y),max(y),length=ncol(mat)),
			z=mat,ticktype="detailed",col=color[facetcol],phi=phi,theta=theta,
			zlab=zlab,xlab=xlab,ylab=ylab,main=main,axes=TRUE,...)

		# return info
	        return(invisible(list(x=x,y=y,xd=xd,yd=yd,mat=mat,col=color[facetcol])))
	}

	# dynamic 3d plot
	if(plot.type=="persp3d"){
		# the color portion of this code is adapted from the persp() help page
		#par(bg="white")
		nrz<-nrow(mat)
		ncz<-ncol(mat)

		# Create a function interpolating colors in the range of specified colors
        	if(color=="heat"){
            		jet.colors<-colorRampPalette(heat.colors(100))
        	}else if(color=="topo"){
			jet.colors <- colorRampPalette(topo.colors(100)) 
        	}else if(color=="cm"){
            		jet.colors<-colorRampPalette(cm.colors(100))
        	}else if(color=="terrain"){
            		jet.colors<-colorRampPalette(terrain.colors(100))
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		jet.colors<-colorRampPalette(gray(seq(0.1,0.9,length=7)))
		}

		# Generate the desired number of colors from this palette
		nbcol<-100
		color<-jet.colors(nbcol)

		# Compute the z-value at the facet centres
		zfacet<-mat[-1,-1]+mat[-1,-ncz]+mat[-nrz,-1]+mat[-nrz,-ncz]

		# Recode facet z-values into color indices
		facetcol<-cut(zfacet,nbcol)
		facetcol=color[facetcol]

		# this portion is from the persp3d() help page
		nx=nrow(mat)
		ny=ncol(mat)
		col <- rbind(1, cbind(matrix(facetcol, nx-1, ny-1), 1))

		# plot
		persp3d(x=seq(min(x),max(x),length=nrow(mat)),y=seq(min(y),max(y),length=ncol(mat)),
			z=mat,col=col,zlab=zlab,main=main,alpha=alpha,smooth=FALSE,lit=lit,
			xlab=xlab,ylab=ylab,...)

		# return info
		return(invisible(list(x=x,y=y,xd=xd,yd=yd,mat=mat,col=pal)))
	}
}
