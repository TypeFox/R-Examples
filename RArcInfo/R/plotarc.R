#Plots all the arcs imported from an ARC file by get.arcdata

#New: T for new plots
plotarc<-function(arc, new=TRUE, index=NULL, xlim, ylim, ...)
{
	if(is.null(index))
		index<-1:length(arc[[2]])

	ll<-length(index)

	#We only need the list of arcs, not the dataframe with the other data
	arc<-arc[[2]][index]

	if(new==TRUE)
	{
		if(missing(xlim) || missing(ylim) )
		{
			#Calculate the boundary

			x<-arc[[1]][[1]]
			y<-arc[[1]][[2]]
			l<-as.integer(length(x))

			xmin<-min(x)
			xmax<-max(x)
			ymin<-min(y)
			ymax<-max(y)

			nxmin<-min(x)
			nxmax<-max(x)
			nymin<-min(y)
			nymax<-max(y)


			for(i in 2:ll)
			{
				x<-arc[[i]][[1]]
				y<-arc[[i]][[2]]
				l<-as.integer(length(x))

				nxmin<-min(x)
				nxmax<-max(x)
				nymin<-min(y)
				nymax<-max(y)

				if(nxmin<xmin) xmin<-nxmin
				if(nymin<ymin) ymin<-nymin
				if(nxmax>xmax) xmax<-nxmax
				if(nymax>ymax) ymax<-nymax
			}
		}
		else #if(!exists("xlim") || !exists("ylim") )
		{
			xmin<-xlim[1]
			xmax<-xlim[2]
			ymin<-ylim[1]
			ymax<-ylim[2]
		}

		range<-max( c(xmax-xmin, ymax-ymin) )/2

		xmean<-(xmin+xmax)/2
		ymean<-(ymin+ymax)/2


		#Set aspect ratio and display plotting window

		par.in <- par(no.readonly = TRUE)
#		on.exit(par(par.in))

		plot.dim<-c(xmax-xmin, ymax-ymin)
		print(min(par.in$pin)
		                * par.in$fin / max(par.in$fin)
				                * (plot.dim) / max(plot.dim))
		par(pin = min(par.in$pin) 
		* par.in$fin / max(par.in$fin)
		* (plot.dim) / max(plot.dim))

		plot(arc[[1]][[1]], arc[[1]][[2]], xlim=c(xmean-range,xmean+range), ylim=c(ymean-range, ymean+range), type="n", ...)

	}##if(new)

	for (i in 1:ll)
	{
		lines(arc[[i]][[1]], arc[[i]][[2]], ...)
	}

}
