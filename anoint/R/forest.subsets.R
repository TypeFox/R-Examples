percent.inner <- 0.9

inner.view <- viewport(width=percent.inner,height=percent.inner,
										layout=grid.layout(nrow=1,ncol=2,widths=1/2))
	
drawgrid <- function(OnMatrix,on.fill="blue",off.fill="white"){
	
	nr <- nrow(OnMatrix)
	nc <- ncol(OnMatrix)
	
	width <- 1/nc
	height <- 1/nr
	
	y <- 1

	for(i in 1:nr){
		xpos <- seq(0,1,by=width)
		xpos <- xpos[-length(xpos)]
		
		colors <-  off.fill
		on <- OnMatrix[i,]
		colors[on] <- on.fill
		grid.rect(x=xpos,y=y,width=width,height=height,gp=gpar(fill=colors),just=c("left","top"))
		y <- y-height
	}

}

# PLACES MINIMUM AT 0 AND MAXIMUM AT 1
convert.npc <- function(y=1:5, limits=NULL, clip.replace=NULL){

	if(is.null(limits))	 limits <- range(y)

	if(is.null(clip.replace))
		y[y<limits[1]|y>limits[2]] <- NA
	else
		y[y<limits[1]|y>limits[2]] <- clip.replace
		
	if(any(is.na(y)))
		warning("Values not within limits will be clipped.")

(y-limits[1])/diff(limits)
}


draweffects <- function(effects, confint=TRUE, lower=NULL, upper=NULL,segments.gpar=NULL,
														xlimits=NULL, reject, fill, legend=TRUE, subgroup=FALSE){
	
	y <- (0:(length(effects)-1)+1:length(effects))/2
	y.npc <- convert.npc(y, c(0,length(effects)))
	y.npc <- sort(y.npc, decreasing=TRUE)
	if(!subgroup)
		x.npc <- convert.npc(c(effects,1), xlimits)
	else
		x.npc <- convert.npc(c(effects,0), xlimits)
	col <- (c("white",fill))[factor(reject)]

	if(confint){
		lower.npc <- convert.npc(lower, xlimits, clip.replace=xlimits[1])
		upper.npc <- convert.npc(upper, xlimits, clip.replace=xlimits[2])
		if(is.null(segments.gpar)) segments.gpar <- get.gpar()
		grid.segments(x0=lower.npc, x1=upper.npc, y0=y.npc, y1=y.npc, gp=segments.gpar)
	}
	grid.circle(x=x.npc[-length(x.npc)], y=y.npc, r=1/max(c(5*length(effects),30)),
					gp=gpar(fill=col))
	
	grid.lines(x=c(0,1), y=c(0,0)) # HORIZONTAL LINE
	grid.lines(x=rep(x.npc[length(x.npc)],2), y=c(0,1)) # CENTER LINE
	
	# LEGEND
	if(legend){
		grid.circle(x=c(0.95, 0.95), y=c(1,0.97), r=1/50,
					gp=gpar(fill=c("white",fill)))
		grid.text(x=0.98, y=1.03, label="Significance")			
	    grid.text(x=c(0.975, 0.975), y=c(1,0.97), label=c("No","Yes"), just=c("left","centre"))
	}
}

effectaxis <- function(xlimits, xat=NULL, percent.inner,effects.text=NULL,effects.axis=NULL){

	if(is.null(xat)) xat <- seq(xlimits[1], xlimits[2], length=6)
	
	x <- convert.npc(xat, xlimits)
	
	grid.segments(x0=x, x1=x, y0=-1/100, y1=1/100)

	popViewport(1)
	
	pushViewport(viewport(width=percent.inner))
	
	if(is.null(effects.text)) effects.text <- get.gpar()	
	if(is.null(effects.axis)) effects.axis <- get.gpar()
	
	grid.text(label=format(xat, ns=2, dig=2), x=x, y=1/40, just=c("centre","top"),gp=effects.axis)
	grid.text(y=1, x=0.5, label="Interaction Effect",gp=effects.text)
}
 

subgroupaxis <- function(labels=NULL, 
										percent.inner, 
										subgroup.title=NULL,
									    subgroup.text=NULL, 
									    subgroup.axis=NULL,
									    xlimits=c(0,2)){

	popViewport(2)
	pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
	pushViewport(viewport(width=percent.inner))
	x <- convert.npc((1:length(labels)+0:(length(labels)-1))/2, c(0, length(labels)))

	if(is.null(subgroup.text)) subgroup.text <- get.gpar()
	if(is.null(subgroup.axis)) subgroup.axis <- get.gpar()
	
	grid.text(y=1/75,x=x,label=labels,rot=30, gp=subgroup.axis)
	grid.text(y=1,x=0.5,subgroup.title,gp=subgroup.text)
}


forest.subsets <- function(object,
										index=1:(min(length(object$interaction),30)),
										labels=NULL,
										exclude.fill="white", 
										include.fill="grey30", 
										signif.fill = "red",
										percent.inner=0.9,
										xlimits=NULL,
										legend=TRUE,
										subgroup.text=NULL,
										subgroup.axis=NULL,
										subgroup.title="Included Covariates",
										effects.text=NULL,
										effects.axis=NULL,
										confint=TRUE,
										segments.gpar=NULL,
										subgroup=FALSE){


	grid.newpage()
	
	inner.view <- viewport(width=percent.inner,height=percent.inner,
										layout=grid.layout(nrow=1,ncol=2,widths=1/2))

	pushViewport(inner.view)

	pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
	
	pushViewport(viewport(width=percent.inner, height=percent.inner))

	drawgrid(object$include.exclude.matrix[index,], on.fill=include.fill, off.fill=exclude.fill)

	popViewport(2)
	pushViewport(viewport(layout.pos.col=2,layout.pos.row=1))
	pushViewport(viewport(width=percent.inner, height=percent.inner))

	if(is.null(xlimits)){
		if(confint)
			xlimits <- range(c(object$lower[index],object$upper[index]))+c(-0.1,0.1)
		else
			xlimits <- range(object$interaction[index])+c(-0.1,0.1)
	}

	draweffects(object$interaction[index], confint=confint,
									lower=object$lower[index], upper=object$upper[index], 																					segments.gpar=segments.gpar,
								    xlimits=xlimits, reject = object$reject[index], 
									fill=signif.fill, legend=legend, subgroup=subgroup)

	effectaxis(xlimits=xlimits, percent.inner=percent.inner,
							effects.text=effects.text, 	effects.axis=effects.axis)

	if(is.null(labels)) labels <- object$covariates
	
	subgroupaxis(labels, percent.inner=percent.inner,subgroup.title=subgroup.title,
							subgroup.text=subgroup.text, 	subgroup.axis=subgroup.axis)

}


