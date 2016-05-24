plotRaw3d.fnc<-function(data=NULL,
			response=NULL,
			pred=NULL,
			intr=NULL,
			xy=TRUE,
			color="topo",
			zlim=NULL,
			xlab=NULL,
			ylab=NULL,
			zlab=NULL,
			main=NULL,
			shift=0,
                        scale=1,
			plot.type="contour",
			add=FALSE,
			alpha=1,
			theta=30,
			phi=30,
			ticktype="detailed",
			contourstepsize=1,
                        legend.args=NULL,
                        ...
){
	if(is.null(data))stop("please specify a data frame\n")
	if(is.null(response))stop("please specify a response variable\n")
	if(is.null(pred))stop("please specify a predictor\n")
	if(is.null(intr))stop("please specify an interacting predictor\n")

	# set labels if NULL
	if(is.null(xlab)){
		xlab=pred
	}

	if(is.null(ylab)){
		ylab=intr
	}

	if(is.null(zlab)){
		zlab=response
	}

	if(is.null(main)){
		if(plot.type=="contour"){
			main=zlab
		}else{
			main=""
		}
	}

	# get average
	x<-tapply(data[,response],list(data[,pred],data[,intr]),function(x)mean(x,na.rm=TRUE))
	x[is.na(x)]<-0

        x<-x*scale+shift
	
	if(is.null(zlim[1])){
		zlim=range(x)
	}

	if(plot.type[1]=="image.plot"){
    	contourlevels = seq(zlim[1], zlim[2], by=contourstepsize)
	                                        
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
	    stop("color scheme not recognised")
	  }

	  image.plot(x=as.numeric(rownames(x)),y=as.numeric(colnames(x)),z=x,col=pal,
	    zlim=zlim,xlab=xlab,ylab=ylab,main=zlab,...)
          contour(x=as.numeric(rownames(x)),y=as.numeric(colnames(x)),z=x,
	    col=con.col,zlim=zlim,add=TRUE,levels=round(contourlevels,2),...)

          return(invisible(list(z=x,col=pal)))
      }else if(plot.type[1]=="persp3d"){

		dev.new()

		# the color portion of this code is adapted from the persp() help page
		#par(bg="white")
		nrz<-nrow(x)
		ncz<-ncol(x)
	
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
	        }else{
			stop("color scheme not recognised")
		}
	
		# Generate the desired number of colors from this palette
		nbcol<-100
		color<-jet.colors(nbcol)
	
		# Compute the z-value at the facet centres
		zfacet<-x[-1,-1]+x[-1,-ncz]+x[-nrz,-1]+x[-nrz,-ncz]
	
		# Recode facet z-values into color indices
		facetcol<-cut(zfacet,nbcol)
		facetcol=color[facetcol]
	
		# this portion is from the persp3d() help page
		nx=nrz
		ny=ncz
		col <- rbind(1, cbind(matrix(facetcol, nx-1, ny-1), 1))


		if(add){
			box=FALSE
			axes=FALSE
			main=""
			xlab=""
			ylab=""
			zlab=""
		}else{
			box=TRUE
			axes=TRUE
		}
	
		if(xy){
			persp3d(x=as.numeric(rownames(x)),y=as.numeric(colnames(x)),z=x,
				xlab=xlab,ylab=ylab,zlab=zlab,main=main,col=col,zlim=zlim,
				smooth=FALSE,add=add,alpha=alpha,box=box,axes=axes)
		}else{
			persp3d(z=x,xlab=xlab,ylab=ylab,zlab=zlab,main=main,col=col,
				zlim=zlim,smooth=FALSE,add=add,alpha=alpha,box=box,axes=axes)
		}

		dev.off()

		return(invisible(list(z=x,col=col)))
	}else if(plot.type[1]=="persp"){
		# the color portion of this code is adapted from the persp() help page
		#par(bg="white")
		nrz<-nrow(x)
		ncz<-ncol(x)
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
        		}else{
			stop("color scheme not recognised")
		}
		# Generate the desired number of colors from this palette
		nbcol<-100
		color<-jet.colors(nbcol)
		# Compute the z-value at the facet centres
		zfacet<-x[-1,-1]+x[-1,-ncz]+x[-nrz,-1]+x[-nrz,-ncz]
		# Recode facet z-values into color indices
		facetcol<-cut(zfacet,nbcol)

		persp(x=as.numeric(rownames(x)),y=as.numeric(colnames(x)),z=x,xlab=xlab,
			ylab=ylab,zlab=zlab,main=main,col=color[facetcol],zlim=zlim,theta=theta,
			phi=phi,ticktype=ticktype)
		return(invisible(list(z=x,col=color[facetcol])))
	}else if(plot.type[1]=="contour"){
			contourlevels = seq(zlim[1], zlim[2], by=contourstepsize)
			
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
				stop("color scheme not recognised")
			}

			image(x=as.numeric(rownames(x)),y=as.numeric(colnames(x)),z=x,col=pal,
				zlim=zlim,xlab=xlab,ylab=ylab,main=zlab)
			contour(x=as.numeric(rownames(x)),y=as.numeric(colnames(x)),z=x,
				col=con.col,zlim=zlim,add=TRUE,levels=round(contourlevels,2))
			box()
			return(invisible(list(z=x,col=pal)))
	}else{
		stop("Plot type unrecognizable!\n")
	}
}
