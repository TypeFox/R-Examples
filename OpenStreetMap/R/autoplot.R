

#' Plots an open street map tile using ggplot2
#' @param data an osmtile
#' @param plot if false only the annotation_raster is returned
#' @param ... not used
#' @method autoplot osmtile
autoplot.osmtile <- function(data,plot=FALSE,...){
	a <- b <- NULL
	x <- data
	p1 <- x$bbox$p1
	p2 <- x$bbox$p2
	yres <- x$yres
	xres <- x$xres
	rast <- as.raster(matrix(x$colorData, nrow=x$xres, 
					byrow = TRUE))
	annot <- ggplot2::annotation_raster(rast,
			p1[1] - .5*abs(x$bbox$p1[1]-x$bbox$p2[1])/yres,
			p2[1] + .5*abs(x$bbox$p1[1]-x$bbox$p2[1])/yres,
			p2[2] - .5*abs(x$bbox$p1[2]-x$bbox$p2[2])/xres,
			p1[2] + .5*abs(x$bbox$p1[2]-x$bbox$p2[2])/xres)
	if(plot)
		ggplot2::ggplot(ggplot2::aes(x=a,y=b),data=data.frame(a=p1[1],b=p1[2])) + 
				annot + ggplot2::expand_limits(x = c(p1[1],p2[1]),
						y=c(p2[2],p1[2]))
	else
		annot
}

#' Plot an open street map using ggplot2
#' @param data an OpenStreetMap object
#' @param expand if true the plotting bounds are expanded to the bounding box
#' @param ... not used
#' @examples \dontrun{
#' require(maps)
#' require(ggplot2)
#' 
#' mp <- openmap(c(53.38332836757155,-130.517578125),
#' 		c(15.792253570362446,-67.939453125),4,'stamen-watercolor')
#' mp_bing <- openmap(c(53.38332836757155,-130.517578125),
#' 		c(15.792253570362446,-67.939453125),4,'bing')
#' states_map <- map_data("state")
#' states_map_merc <- as.data.frame(
#' 		projectMercator(states_map$lat,states_map$long))
#' states_map_merc$region <- states_map$region
#' states_map_merc$group <- states_map$group
#' crimes <- data.frame(state = tolower(rownames(USArrests)), USArrests)
#' 
#' p <- autoplot(mp,expand=FALSE) + geom_polygon(aes(x=x,y=y,group=group),
#' 		data=states_map_merc,fill="black",colour="black",alpha=.1) + theme_bw()
#' print(p)
#' p <- autoplot(mp_bing) + geom_map(aes(x=-10000000,y=4000000,map_id=state,fill=Murder),
#' 		data=crimes,map=states_map_merc)
#' print(p)
#' }
#' @method autoplot OpenStreetMap
autoplot.OpenStreetMap <- function(data, expand=TRUE, ...){
	x <- y <- NULL
	x <- data
	p1 <- x$bbox$p1
	p2 <- x$bbox$p2
	p <- ggplot2::ggplot(ggplot2::aes(x=x,y=y),data=data.frame(x=(p1[1]+p2[1])/2,
					y=(p1[2]+p2[2])/2))
	if(expand)
		p <- p + ggplot2::expand_limits(x = c(p1[1],p2[1]),
						y=c(p2[2],p1[2])) + 
				ggplot2::scale_x_continuous(expand=c(0,0)) + 
				ggplot2::scale_y_continuous(expand=c(0,0))
	for(tile in x$tiles){
		p <- p + autoplot.osmtile(tile)
	}
	p + ggplot2::coord_equal()
}


