plotpoly <-function(arc,bnd,pal,index=NULL,col, border=NULL, xratio=1, yratio=1, ...)
{
	if(is.null(index))
		index<-1:length(pal[[1]][[1]])

	lindex<-length(index)
	col<-rep(col,length.out=lindex)	
	
	if(is.null(border))
		border<-rep("black", length(col))
	else
		border<-rep(border, length.out=lindex)


#Set aspect ratio and display plotting window

	par.in <- par(no.readonly = TRUE)
#	on.exit(par(par.in))

	plot.dim<-c(bnd[3]-bnd[1], bnd[4]-bnd[2])

	rdib<-min(par.in$pin[1]/plot.dim[1],par.in$pin[2]/plot.dim[2])
	plotreg<-c(bnd[[1]],bnd[[1]]+plot.dim[[1]]*rdib,
			bnd[[2]],bnd[[2]]+plot.dim[[2]]*rdib)
	par(pin=c(xratio,yratio)*plot.dim*rdib, usr=plotreg)

	plot((bnd[1]+bnd[3])/2,(bnd[2]+bnd[4])/2,xlim=c(bnd[1],bnd[3]),
			ylim=c(bnd[2],bnd[4]),type="n", ...)	


#First, we just select the arcs we will need. This will speed up
#this function
	palindex<-match(index, pal[[1]]$PolygonId)
	arcindex<-as.vector(sapply(pal[[2]][palindex],function(X){X[[1]]}))
	arcindex<-sort(unique(abs(unlist(arcindex))))[-1]
	
	arc1<-list(arc[[1]][1:7][arcindex,],arc[[2]][arcindex])
	pal1<-list(pal[[1]][1:6][palindex,],pal[[2]][palindex])


	for(i in 1:lindex)
	{

#Now, we plot the polygons

		#This gives us the list of arcs 
		p<-match(index[i], pal1[[1]]$PolygonId)
		p<-pal1[[2]][[p]][[1]]


		#And here we get the "rows" of them in the variable arc
		absp<-abs(p)
		a<-match(absp, arc1[[1]]$ArcId)

		#When p[[j]]==0 it means that the previous arcs build a closed polygon
		#and it can be plotted

		#Default values for x and y coordinates.
		x<-c()
		y<-c()
	
		
		for( j in 1:length(p) )
		{	
	
			if(p[[j]]>0)
                	{
        	                x<-c(x, arc1[[2]][[ a[j] ]][[1]])
       	        	        y<-c(y, arc1[[2]][[ a[j] ]][[2]])
        	        }
       	        	else if(p[[j]]<0)
                	{
				x<-c(x, rev(arc1[[2]][[ a[j] ]][[1]]))
				y<-c(y, rev(arc1[[2]][[ a[j] ]][[2]]))
       	        	}
			else
			{
				if(length(x)>0 && length(y)>0)
				{
					if( (x[1]!=x[length(x)]) || (y[1]!=y[length(y)]) )
					{
						print("The polygon isn't closed")
						x<-c(x,x[1])
						y<-c(y,y[1])
						
						idx<-match(index[i], index)
						polygon(x,y,col=col[idx], border=border[idx])
					}
					else
					{
						idx<-match(index[i], index)
						polygon(x,y,col=col[idx], border=border[idx])
					}

					x<-c()
					y<-c()
				}
			}
		}


		if(length(x)>0 && length(y)>0)
		{
			if(x[1]!=x[length(x)] || y[1]!=y[length(y)] )
			{
				print("The polygon isn't closed")
				print(c(j,p[[j]]))
				x<-c(x,x[1])
				y<-c(y,y[1])
			}

			idx<-match(index[i], index)
			polygon(x,y,col=col[idx], border=border[idx])
		}
	}
}
