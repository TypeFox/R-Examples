topoplot.palette<-function(cols, pos = c(0.5,0.5), p.width=0.2, p.height=0.8, horizontal = FALSE, rev= FALSE, palette.lwd=1, palette.bord="black", palette.lim = c(-5,5), labels=TRUE, lab.dist = 1, lab.cex = 1, lab.font = 1, lab.family = "")
{
	# reverse vector if required
	if (rev == T){
		cols=rev(cols)
		palette.lim=rev(palette.lim)
	}
	
	# some adjustments to palette.lim
	palette.lim[palette.lim>0]=paste("+", palette.lim[palette.lim>0], sep="")
	
	#find coordinates
	xlim=par("usr")[1:2]
	ylim=par("usr")[3:4]
	
	tot.x=as.numeric(dist(xlim))
	tot.y=as.numeric(dist(ylim))
	
	#calculate the position of the center of the rectangle (x.axis)
	pos.x=tot.x*pos[1]+xlim[1]
	
	#calculate the position of the center of the rectangle (y.axis)
	pos.y=tot.y*pos[2]+ylim[1]
	
	### find rectangle dimensions
	
	rect.length=p.width*tot.x
	rect.height=p.height*tot.y
	
	x0rect=pos.x-(rect.length/2)	
	x1rect=pos.x+(rect.length/2)
	
	y0rect=pos.y-(rect.height/2)
	y1rect=pos.y+(rect.height/2)
	
	# draw palette horizontally
	
	if (horizontal==TRUE){
	
	#find palette steps length
	p.steps.length=rect.length/length(cols)	

	for (i in 0:(length(cols)-1)){
		
		p.xcoord=c(x0rect+i*p.steps.length, x0rect+i*p.steps.length+p.steps.length, x0rect+i*p.steps.length+p.steps.length, x0rect+i*p.steps.length)		
		p.ycoord=c(y0rect, y0rect, y1rect, y1rect)
		polygon(p.xcoord,p.ycoord, col=cols[i+1], border=cols[i+1])
		
	}
		
	#draw rectangle
	polygon(c(x0rect, x1rect, x1rect, x0rect), c(y0rect, y0rect, y1rect, y1rect), border = palette.bord, lwd = palette.lwd)
	
	#### draw labels when plot is horizontal

	if (labels == TRUE){
	
	# calculate the distance between the label and the rectangle. 
	lab.x.dist=tot.x*lab.dist/100
		
	
	text((x0rect - lab.x.dist) ,  pos.y, labels= bquote(paste(.(palette.lim[1]), " ",  mu, "V")), pos = 2, cex = lab.cex, font = lab.font, family = lab.family) 
	# notice how I retrieve the value of a variable with .()
	text((x1rect + lab.x.dist) ,  pos.y, labels= bquote(paste(.(palette.lim[2]), " ",  mu, "V")), pos = 4, cex = lab.cex, font = lab.font, family = lab.family)

	
	}

	
	}
	
	
	# draw palette vertically
	
	if (horizontal==FALSE){
		#find palette steps length
	p.steps.length=rect.height/length(cols)	

	for (i in 0:(length(cols)-1)){
		
		p.ycoord=c(y0rect+i*p.steps.length, y0rect+i*p.steps.length+p.steps.length, y0rect+i*p.steps.length+p.steps.length, y0rect+i*p.steps.length)		
		p.xcoord=c(x0rect, x0rect, x1rect, x1rect)
		polygon(p.xcoord,p.ycoord, col=cols[i+1], border=cols[i+1])
		
	}
		
	#draw rectangle when it is vertical
	polygon(c(x0rect, x1rect, x1rect, x0rect), c(y0rect, y0rect, y1rect, y1rect), border = palette.bord, lwd = palette.lwd)
	
	#### draw labels when plot is vertical
	if (labels == TRUE){
	
	# calculate the distance between the label and the rectangle. 
	lab.y.dist=tot.y*lab.dist/50
		
	
	text(pos.x , y0rect-lab.y.dist, labels= bquote(paste(.(palette.lim[1]), " ",  mu, "V")),  cex = lab.cex, font = lab.font, family = lab.family) 
	# notice how I retrieve the value of a variable with .()
	text(pos.x, y1rect+lab.y.dist, labels= bquote(paste(.(palette.lim[2]), " ",  mu, "V")),  cex = lab.cex, font = lab.font,  family = lab.family)
	
	invisible(lab.y.dist)
	
	}

			
	}
	
	
	
	
}