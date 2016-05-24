

spatialColoredPoints<-function(x, color_var, pch=1, palette , legend.loc="bottomleft",
		legend.title,...){
	if(length(color_var)==1 && is.character(color_var)){
		if(missing(legend.title))
			legend.title = color_var
		color_var <- slot(x,"data")[,color_var]
	}else if(missing(legend.title))
		legend.title <- deparse(substitute(color_var))
	if(is.character(color_var))
		color_var <- as.factor(color_var)
	if(is.factor(color_var)){
		if(missing(palette))
			palette <- function(n)rainbow(n,v=.8)
		nlevs <- length(levels(color_var))
		if(length(palette(3))<=1){
			pal <- palette(1:nlevs/nlevs)
			palette <- manual_pal(pal)
		}
		clrs <- palette(nlevs)
		org <- levels(color_var) 
		levels(color_var) <- clrs
		color_var <- as.character(color_var)
		plot(x,col=color_var,add=TRUE,pch=pch,...)
		legend(legend.loc,,org,col=clrs,pch=pch,title=legend.title)
	}else{
		#color_var <- as.numeric(color_var)
		if(missing(palette))
			palette <- gradient_n_pal(heat.colors(20))
		if(length(palette(3))>1){
			pal <- palette(9)
			palette <- gradient_n_pal(pal)
		}
		clrs <- cscale(color_var,palette)

		qnt <- rev(quantile(color_var,type=1,na.rm=TRUE))
		qntclrs <- sapply(qnt,function(x)clrs[which(color_var==x)[1]])
		#repl <- clrs[cv]
		#leg.col <- clrs[c(1,25,50,100)]
		#leg.val <- c(.01,.25,.50,1)*(max(color_var,na.rm=TRUE)-min(color_var,na.rm=TRUE)) + min(color_var,na.rm=TRUE)
		#leg.val <- format(leg.val,digits=3)
		plot(x,col=clrs,add=TRUE,pch=pch,...)
		legend(legend.loc,,formatC(qnt,digits=4),col=qntclrs,
				pch=pch,title=legend.title, bg = "white")
	}
}




spatialBubblePlot <- function(x,z,minRadius=.01,
		maxRadius=.05,color="#F75252", ...){
	if(length(z)==1 && is.character(z))
		z <- slot(x,"data")[,z]
	mat <- coordinates(x)
	z <- (z-min(z, na.rm=TRUE)) 	
	z <- sqrt(z/pi)
	z <- z / max(z, na.rm=TRUE)
	r <- z*(maxRadius-minRadius) + minRadius;
	dd <- data.frame(x=mat[,1],y=mat[,2],r=r);
	dd <- dd[order(-r),];
	symbols(dd$x, dd$y, circles=dd$r,
			inches=maxRadius,add=TRUE, fg="white", bg=color, ...);
}


spatialTextPlot <- function(x,text,...){
	coord <- coordinates(x)
	if(length(text)==1 && is.character(text))
		text <- slot(x,"data")[text]
	text(coord[,1],coord[,2],text,...)
}



spatialChoropleth <- function (sp, color ,quantileBin=FALSE, palette, alpha=1,
		main = NULL, sub = "", legend.loc = "bottomleft", 
		legend.title ,add=TRUE,border="transparent", ...) 
{
	cname <- deparse(substitute(color))
	if(is.character(color) && length(color)==1){
		cname <- color
		color <- slot(sp,"data")[,color]
	}
	if(is.numeric(color) && quantileBin>0){
		if(quantileBin==1)
			quantileBin <- 5
		color <- cut2(color,g=quantileBin)
		if(missing(palette))
			palette <- heat.colors
	}
	if(is.character(color))
		color <- as.factor(color)
	if(is.factor(color)){
		if(missing(palette))
			palette <- function(n)rainbow(n,v=.8)
		if(length(palette(3))<=1){
			s <- palette(seq(from=0,to=1,length.out=length(levels(color))))
			palette <- manual_pal(s)
		}
		clrs <- palette(length(levels(color)))
		org <- levels(color) 
		levels(color) <- clrs
		color <- as.character(color)
		legLabel <- rev(org)
		legColors <- rev(clrs)
	}else{
		#color_var <- as.numeric(color_var)
		if(missing(palette))
			palette <- gradient_n_pal(heat.colors(20))
		if(length(palette(3))>1)
			palette <- gradient_n_pal(palette(9))
		clrs <- cscale(color,palette)
		qnt <- rev(quantile(color,type=1,na.rm=TRUE))
		qntclrs <- sapply(qnt,function(x)clrs[which(color==x)[1]])
		legLabel <- formatC(qnt,digits=4)
		legColors <- qntclrs
		color <- clrs
	}
	if(missing(legend.title)){
		legend.title <- cname
	}
	plot(sp,col=alpha(color,alpha),add=add,border=border,...)
	title(main = main, sub = sub)
	legend(legend.loc, legend = legLabel, fill = legColors, 
			bty = "o", title = legend.title, bg = "white")
}

