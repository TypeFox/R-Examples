##' spplot1 function
##'
##' A function to provide spplot-like plotting capability but NOT using trellis graphics. This function also acts as an interface for fast
##' plotting of SpatialPolygonsDataFrame or SpatialPixelsDataFrame objects using leaflet HTML plotting capabilities to get zoomable plots
##' with real-world context: transformation to the correct projection is done automatically.
##'
##' See \url{http://leaflet-extras.github.io/leaflet-providers/preview/} for examples of  leaflet templates.
##'
##' Instructions on installing the leaflet R package are available from  \url{https://rstudio.github.io/leaflet/}
##'
##' @param x a SpatialPolgonsDataFrame or a SpatialPointsDataFrame 
##' @param what the name of the variable to plot
##' @param palette the palette, can either be a vector of names of colours, or a vector of colours produced for example by the brewer.pal function.
##' @param breaks optional breaks for the legend, a vector of length 1 + length(palette)
##' @param legpos the position of the legend, options are 'topleft', 'topright', 'bottomleft', 'bottomright'
##' @param fun an optional function of the data to plot, default is the identity function
##' @param include.lowest see ?cut
##' @param bty see ?legend
##' @param bg see ?legend
##' @param printlegend logical: print the legend?
##' @param bw Logical. Plot in black/white/greyscale? Default is to produce a colour plot. Useful for producing plots for journals that do not accept colour plots.
##' @param useLeaflet whether to use leaflet to produce a zoomable map this requires the leaflet package, available by issuing the command "devtools::install_github('rstudio/leaflet')"
##' @param urltemplate template for leaflet map background, default is urlTemplate('Stamen-Toner'), but any valid web address for leaflet templates will work here. See ?urlTemplate.
##' @param fillOpacity see ?addPolygons
##' @param legendOpacity see opacity argument in function addLegend
##' @param OSMbg optional OpenStreetMap background to add to plot, obtain this using the function getBackground
##' @param leafletLegend logical, display the leaflet legend?
##' @param ... other arguments to be passed to plot
##' @return either produces a plot or if useLeaflet is TRUE, returns a leaflet map widget to which further layers can be added
##' @seealso \link{urlTemplate}, \link{getBackground}, brewer.pal
##' @export

spplot1 <- function(x,
					what,
					palette=brewer.pal(5,"Oranges"),
					breaks=NULL,
					legpos="topleft",
					fun=identity,
					include.lowest=TRUE,
					bty="n",
					bg=NULL,
					printlegend=TRUE,
					bw=FALSE,
					useLeaflet=FALSE,
					urltemplate=urlTemplate("Stamen_Toner"),
					fillOpacity = 0.5,
					legendOpacity = 0.5,
					OSMbg=NULL,
					leafletLegend=TRUE,
					...){

	if(bw){
		palette <- brewer.pal(5,"Greys")	
	}
	n <- length(palette)
	if(!is.null(breaks)){
		if(length(breaks)!=(n+1)){
			stop("Breaks must be a vector of length 1 + length(palette).")
		}
	}

	if(is.numeric(x[[what]])){
		tp <- fun(x[[what]]) # tp = to plot
		if(is.null(breaks)){
			cutt <- cut(tp,n,include.lowest=include.lowest)
		}
		else{
			cutt <- cut(tp,breaks,include.lowest=include.lowest)	
		}
		lvls <- levels(cutt)
		cols <- palette[as.numeric(cutt)]
	}
	else if(is.factor(x[[what]])){
		cutt <- as.numeric(x[[what]])
		lvls <- levels(x[[what]])
		if(length(lvls)!=length(palette)){
			stop(paste("Number of levels in factor variable, ",what,", must have the same number of colours in palette.",sep=""))
		}
		cols <- palette[cutt]
	}
	else{
		stop("Plotting variable must either be of class numeric, integer or factor.")
	}

	
	
	if(!useLeaflet){
		if(is.null(OSMbg)){
			sp::plot(x,col=cols,...)
		}
		else{
			plot(OSMbg,...)
			cols <- adjustcolor(cols,alpha.f=0.5)
			palette <- adjustcolor(palette,alpha.f=0.5)
			x <- spTransform(x,CRS("+init=epsg:3857"))
			sp::plot(x,col=cols,add=TRUE,...)
		}

		if(printlegend){
			legend(legpos,pch=rep(15,n),col=palette,legend=lvls,bty=bty,bg=bg)
		}
	}
	else{
		m <- NULL
		
		s <- "require(leaflet)
		{
			if(inherits(x,'SpatialPolygonsDataFrame')){
				ns <- spTransform(x,CRS('+init=epsg:4326'))
				m <- leaflet() %>%
				addTiles(urlTemplate=urltemplate) %>%
				addPolygons(data=ns,color=cols,fillColor=cols,weight=0,fillOpacity = fillOpacity,stroke=FALSE) %>%
				addLegend(position='topright', labels = lvls, colors = palette,opacity=legendOpacity)
			}
			else if(inherits(x,'SpatialPointsDataFrame')){
				ns <- spTransform(x,CRS('+init=epsg:4326'))
				m <- leaflet() %>%
				addTiles(urlTemplate=urltemplate) %>%
				addCircleMarkers(data=ns,color=cols,stroke=FALSE) %>%
				addLegend(position='topright', labels = lvls, colors = palette,opacity=legendOpacity)
			}
			else{
				stop('Leaflet mapping for this kind of object not supported at present')
			}
		}
		cat('Leaflet map returned, type name of object to visualise in a web browser.\n')"

		if(!leafletLegend){
			s <- "require(leaflet)
			{
				if(inherits(x,'SpatialPolygonsDataFrame')){
					ns <- spTransform(x,CRS('+init=epsg:4326'))
					m <- leaflet() %>%
					addTiles(urlTemplate=urltemplate) %>%
					addPolygons(data=ns,color=cols,fillColor=cols,weight=0,fillOpacity = fillOpacity,stroke=FALSE) 
				}
				else if(inherits(x,'SpatialPointsDataFrame')){
					ns <- spTransform(x,CRS('+init=epsg:4326'))
					m <- leaflet() %>%
					addTiles(urlTemplate=urltemplate) %>%
					addCircleMarkers(data=ns,color=cols,stroke=FALSE) 
				}
				else{
					stop('Leaflet mapping for this kind of object not supported at present')
				}
			}
			cat('Leaflet map returned, type name of object to visualise in a web browser.\n')"
		}

		eval(parse(text=s))

		return(m)

		
	}
}

##' spplot_compare function
##'
##' A function to compare two SpatialPolgonsDataFrame or SpatialPointsDataFrame objects using a unified legend for the variable 
##' of interest in both
##'
##' @param x a SpatialPolgonsDataFrame or a SpatialPointsDataFrame 
##' @param y a SpatialPolgonsDataFrame or a SpatialPointsDataFrame 
##' @param what the name of the variable from x to plot 
##' @param what1 the name of the variable from y to plot. default is to plot the variable of the same name
##' @param palette the palette, can either be a vector of names of colours, or a vector of colours produced for example by the brewer.pal function. 
##' @param legpos the position of the legend, options are 'topleft', 'topright', 'bottomleft', 'bottomright' 
##' @param border see ?spplot 
##' @param fun an optional function of the data to plot, default is the identity function 
##' @param t1 title for the plot of x 
##' @param t2 title for the plot of y
##' @param bw Logical. Plot in black/white/greyscale? Default is to produce a colour plot. Useful for producing plots for journals that do not accept colour plots.
##' @param ... other arguments to be passed to the plot function
##' @return produces a plot comparing x[[what]] and y[[what1]]
##' @export

spplot_compare <- function(x,y,what,what1=what,palette=brewer.pal(9,"Oranges"),legpos="topleft",border=NA,fun=identity,t1="",t2="",bw=FALSE,...){
	if(bw){
		palette <- brewer.pal(5,"Greys")	
	}
	n <- length(palette)
	tpx <- fun(x[[what]])
	tpy <- fun(y[[what1]])
	tp <- c(tpx,tpy)
	cutt <- cut(tp,n)
	lvls <- levels(cutt)
	cols <- palette[as.numeric(cutt)]
	colsx <- cols[1:length(tpx)]
	colsy <- cols[(length(tpx)+1):length(tp)]

	#dev.new(height=6,width=12)
	par(mfrow=c(1,2))

	sp::plot(x,col=colsx,border=border,...)
	legend(legpos,pch=rep(15,n),col=palette,legend=lvls,bty="n",cex=0.75)
	title(t1)

	sp::plot(y,col=colsy,border=border,...)
	#legend(legpos,pch=rep(15,n),col=palette,legend=lvls,bty="n",cex=0.75)
	title(t2)

	par(mfrow=c(1,1))
}


##' getBackground function
##'
##' A function to 
##'
##' @param poly X 
##' @param type X 
##' @return ...
##' @export

getBackground <- function(poly,type="stamen-toner"){
	poly <- spTransform(poly,CRS("+init=epsg:4326"))
	bb <- bbox(poly)
	map <- openmap(upperLeft=c(lat=bb[2,2],lon=bb[1,1]), lowerRight=c(lat=bb[2,1],lon=bb[1,2]),type="stamen-toner")
	return(map)
}

