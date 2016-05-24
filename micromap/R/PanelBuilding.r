labels_build <- function(pl, p, DF, att){ 

	DF$tmp.labels <- DF[,unlist(att[[p]]$panel.data)]

	 #** The defaults are set up so that the text will, in theory, cover the entire 
	 #** width of the panel, whiel maximizing text height as well
	 #** With the way the defaults are set up. default text size needs to be doubled
	tmp.tsize <- att[[p]]$text.size*4
	DF$tmp.y <- -DF$pGrpOrd*att[[p]]$text.size
	tmp.limsy <- c((min(-DF$pGrpOrd)-.5)*att[[p]]$text.size, 
			(max(-DF$pGrpOrd)+.5)*att[[p]]$text.size)

	mln <- max(nchar(as.character(DF$tmp.labels)))
	  if (att[[p]]$align=='right') {
			tmp.limsx <- c(-mln,0)
			DF$tmp.adj <- 1
	    }
	  if (att[[p]]$align=='left') {
			tmp.limsx <- c(0,mln)
			DF$tmp.adj <- 0
	    }
	  if (att[[p]]$align=='center') {
			tmp.limsx <- c(-mln/2,mln/2)
			DF$tmp.adj <- .5
	    }


	  #################################
	  #################################
	  #*** insert space for median row	
	  tmp.median.limsy <- NULL

	  if(att$median.row) {
	    if(!any(DF$median)) DF <- rbind(DF, transform(DF[1,], pGrpOrd=1, pGrp=att$m.pGrp, median=TRUE, 
							rank='', tmp.y=-att[[p]]$text.size,
							tmp.labels=''))


	    tmp.median.limsy <- att[[p]]$text.size*c(-1.5,-.5)
	  }


	  #################################
	  #################################

		
	pl <- 
		ggplot(DF) +
	     	geom_text(aes(x=0, y=tmp.y, label=tmp.labels, 
				hjust=tmp.adj, vjust=.4), 
				family=att[[p]]$text.font, fontface=att[[p]]$text.face, size=tmp.tsize) +
	     	facet_grid(pGrp~., scales="free_y", space="free") +
    		scale_colour_manual(values=att$colors) 

	pl <- plot_opts(p, pl, att)		
	pl <- graph_opts(p, pl, att)		
	pl <- axis_opts(p, pl, att, limsx=tmp.limsx, limsy=c(tmp.limsy,tmp.median.limsy), border=FALSE)

	pl


}


ranks_build <- function(pl, p, DF, att){ 
	tmp.tsize <- att[[p]]$size*4
	DF$tmp.y <- -DF$pGrpOrd*att[[p]]$size
	tmp.limsy <- c((min(-DF$pGrpOrd)-.5)*att[[p]]$size, 
			(max(-DF$pGrpOrd)+.5)*att[[p]]$size)

	mln <- max(nchar(as.character(DF$rank)))
	  if (att[[p]]$align=='right') {
			tmp.limsx <- c(-mln,0)
			DF$tmp.adj <- .9
	    }
	  if (att[[p]]$align=='left') {
			tmp.limsx <- c(0,mln)
			DF$tmp.adj <- .1
	    }
	  if (att[[p]]$align=='center') {
			tmp.limsx <- c(-mln/2,mln/2)
			DF$tmp.adj <- .5
	    }

	  #################################
	  #################################
	  #*** insert space for median row
	  tmp.median.limsy <- NULL

	  if(att$median.row) {
	    if(!any(DF$median)) DF <- rbind(DF, transform(DF[1,], pGrpOrd=1, pGrp=att$m.pGrp, median=TRUE, 
							rank='', tmp.y=-att[[p]]$text.size))


	    tmp.median.limsy <- att[[p]]$text.size*c(-1.5,-.5)
	  }
	
	  #################################
	  #################################


	pl <- 
	 	ggplot(DF) +
		geom_text(aes(x=0, y=tmp.y, label=rank, hjust=tmp.adj, vjust=.4), 
				font=att[[p]]$font, face=att[[p]]$face, size=tmp.tsize) +
     		facet_grid(pGrp~., scales="free_y", space="free") 

	pl <- plot_opts(p, pl, att)		
	pl <- graph_opts(p, pl, att)		
	pl <- axis_opts(p, pl, att, limsx=tmp.limsx, limsy=c(tmp.limsy, tmp.median.limsy), border=FALSE)

	pl
}



dot_legend_build <- function(pl, p, DF, att){ 
	DF$tmp.data <- rep(0, nrow(DF))

	tmp.limsy <- -(range(DF$pGrpOrd) + c(-1,1) * .5)
	tmp.limsx <- c(-.5,.5)

	ncolors <- length(att$colors)

	  #################################
	  #################################
	  #*** insert space for median row	
	  tmp.median.limsy <- NULL

	  if(att$median.row) {
	    if(!any(DF$median)) DF <- rbind(DF, transform(DF[1,], pGrpOrd=1, pGrp=att$m.pGrp, median=TRUE, 
							rank='', tmp.data=0))

	    DF$color[DF$median] <- ncolors
	    tmp.median.limsy <- c(-.5, -1.5)
	  }

	  #################################
	  #################################


	pl <- ggplot(DF) 

	if(att[[p]]$point.border) pl <- pl + geom_point(aes(x=tmp.data, y=-pGrpOrd), 
									colour='black',
									size=att[[p]]$point.size*2.5, 
									pch=att[[p]]$point.type)

	pl <- pl + 
		geom_point(aes(x=tmp.data, y=-pGrpOrd, colour=factor(color)), 
				size=att[[p]]$point.size*2, pch=att[[p]]$point.type) +
	     	facet_grid(pGrp~., scales="free_y", space="free") +
		scale_colour_manual(values=att$colors, guide='none') 

	pl <- plot_opts(p, pl, att)		
	pl <- graph_opts(p, pl, att)		
	pl <- axis_opts(p, pl, att, limsx=tmp.limsx, limsy=c(tmp.limsy,tmp.median.limsy))

	pl
}



dot_build <- function(pl, p, DF, att){ 
	DF$tmp.data <- DF[,unlist(att[[p]]$panel.data)]
	DF$tmp.data1 <- DF[,unlist(att[[p]]$panel.data)]
	DF$tmp.data2 <- c(0, DF$tmp.data1[-nrow(DF)])
	DF$tmp.data3 <- DF$pGrpOrd - 1

	tmp.limsx <- range(DF[,unlist(att[[p]]$panel.data)], na.rm=T)
	if(any(!is.na(att[[p]]$xaxis.ticks))) tmp.limsx <- range(c(tmp.limsx, att[[p]]$xaxis.ticks))
	tmp.limsx <- tmp.limsx + c(-1,1) * diff(tmp.limsx)*.05
	tmp.limsy <- -(range(DF$pGrpOrd) + c(-1,1) * .5)


	if(diff(tmp.limsx)==0) tmp.limsx <- tmp.limsx + c(-.5,.5)

	ncolors <- length(att$colors)
	  #################################
	  #################################
	  #*** insert space for median row	
	  tmp.median.limsy <- NULL

	  if(att$median.row) {
	    if(!any(DF$median)) DF <- rbind(DF, transform(DF[1,], pGrpOrd=1, pGrp=att$m.pGrp, median=TRUE, 
							rank='', tmp.data=median(DF$tmp.data)))

	    DF$color[DF$median] <- ncolors
	    tmp.median.limsy <- c(-.5, -1.5)
	  }

	  #################################
	  #################################


	pl <- ggplot(DF) 

	if(!any(is.na(att[[p]]$add.line))){
		if(length(att[[p]]$add.line.col)==1) att[[p]]$add.line.col <- rep(att[[p]]$add.line.col[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.typ)==1) att[[p]]$add.line.typ <- rep(att[[p]]$add.line.typ[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.size)==1) att[[p]]$add.line.size <- rep(att[[p]]$add.line.size[1], length(att[[p]]$add.line))

		for(j in 1:length(att[[p]]$add.line)) pl <- pl + geom_vline(xintercept = att[[p]]$add.line[j], 
										data=DF, colour=att[[p]]$add.line.col[j], 
										linetype=att[[p]]$add.line.typ[j],
										size=att[[p]]$add.line.size[j])
	  }

	if(att[[p]]$median.line) pl <- pl + geom_vline(aes(xintercept = median(tmp.data)), data = DF, 
							colour = att[[p]]$median.line.col, 
							linetype = att[[p]]$median.line.typ, 
							size = att[[p]]$median.line.size)	

	if(att[[p]]$connected.dots) pl <- pl + geom_segment(aes(x = tmp.data1, y = -pGrpOrd,
								xend = tmp.data2, yend = -tmp.data3),
							data = subset(DF, pGrpOrd>1), 
			            			colour = att[[p]]$connected.col,
							size = att[[p]]$connected.size, 
					            	linetype = att[[p]]$connected.typ)   

	if(att[[p]]$point.border) pl <- pl + geom_point(aes(x=tmp.data, y=-pGrpOrd), 
									colour='black',
									size=att[[p]]$point.size*2.5, 
									pch=att[[p]]$point.type)

	pl <- pl + 
		geom_point(aes(x=tmp.data, y=-pGrpOrd, colour=factor(color)), 
				size=att[[p]]$point.size*2, pch=att[[p]]$point.type) +
	     	facet_grid(pGrp~., scales="free_y", space="free") +
		scale_colour_manual(values=att$colors, guide='none') 

	pl <- plot_opts(p, pl, att)		
	pl <- graph_opts(p, pl, att)		
	pl <- axis_opts(p, pl, att, limsx=tmp.limsx, limsy=c(tmp.limsy,tmp.median.limsy))

	pl
}




dot_cl_build <- function(pl, p, DF, att){ 
	DF$tmp.data1 <- DF[,att[[p]]$panel.data[[1]]]
	DF$tmp.data2 <- DF[,att[[p]]$panel.data[[2]]]
	DF$tmp.data3 <- DF[,att[[p]]$panel.data[[3]]]

	DF$tmp.data4 <- c(0, DF$tmp.data1[-nrow(DF)])
	DF$tmp.data5 <- DF$pGrpOrd - 1

	tmp.limsx <- range(DF[,unlist(att[[p]]$panel.data)], na.rm=T)
	if(any(!is.na(att[[p]]$xaxis.ticks))) tmp.limsx <- range(c(tmp.limsx, att[[p]]$xaxis.ticks))
	tmp.limsx <- tmp.limsx + c(-1,1) * diff(tmp.limsx)*.05
	tmp.limsy <- -(range(DF$pGrpOrd) + c(-1,1) * .5)


	ncolors <- length(att$colors)
	  #################################
	  #################################
	  #*** insert space for median row	
	  tmp.median.limsy <- NULL

	  if(att$median.row) {
	    if(!any(DF$median)) DF <- rbind(DF, transform(DF[1,], pGrpOrd=ncolors+1, pGrp=att$m.pGrp, median=TRUE, rank='', 
								tmp.data1=median(DF$tmp.data1),
								tmp.data2=median(DF$tmp.data2),
								tmp.data3=median(DF$tmp.data3)))
	
	 
	    DF$color[DF$median] <- ncolors
	    tmp.median.limsy <- c(-.5, -1.5)
	  }

	  #################################
	  #################################

	pl <- ggplot(DF) +
		geom_segment(aes(x=tmp.data2, y=-pGrpOrd, xend=tmp.data3, yend=-pGrpOrd, 
				colour=factor(color)), size=att[[p]]$line.width) 


	if(att[[p]]$median.line) pl <- pl + geom_vline(aes(xintercept = median(tmp.data)), data=DF, 
							colour=att[[p]]$median.line.col, 
							linetype = att[[p]]$median.line.typ, 
							size = att[[p]]$median.line.size)

	if(!any(is.na(att[[p]]$add.line))){
		if(length(att[[p]]$add.line.col)==1) att[[p]]$add.line.col <- rep(att[[p]]$add.line.col[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.typ)==1) att[[p]]$add.line.typ <- rep(att[[p]]$add.line.typ[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.size)==1) att[[p]]$add.line.size <- rep(att[[p]]$add.line.size[1], length(att[[p]]$add.line))

		for(j in 1:length(att[[p]]$add.line)) pl <- pl + geom_vline(xintercept = att[[p]]$add.line[j], 
										data=DF, colour=att[[p]]$add.line.col[j], 
										linetype=att[[p]]$add.line.typ[j],
										size=att[[p]]$add.line.size[j])
	  }


	if(att[[p]]$connected.dots) pl <- pl + geom_segment(aes(x = tmp.data1, y = -pGrpOrd,
								xend = tmp.data4, yend = -tmp.data5),
							data = subset(DF, pGrpOrd>1), 
			            			colour = att[[p]]$connected.col,
							size = att[[p]]$connected.size, 
					            	linetype = att[[p]]$connected.typ)   
   
  	if(att[[p]]$point.border) pl <- pl + geom_point(aes(x=tmp.data1, y=-pGrpOrd), 
									size=att[[p]]$point.size*2.5, 
									pch=att[[p]]$point.type,
									data=DF)

	pl <- pl + 
		geom_point(aes(x=tmp.data1, y=-pGrpOrd, colour=factor(color)), 
				size=att[[p]]$point.size*2, pch=att[[p]]$point.type) +
	     	facet_grid(pGrp~., scales="free_y", space="free") +
		scale_colour_manual(values=att$colors, guide='none') 

	pl <- plot_opts(p, pl, att)		
	pl <- graph_opts(p, pl, att)		
	pl <- axis_opts(p, pl, att, limsx=tmp.limsx, limsy=c(tmp.limsy,tmp.median.limsy))
	
	pl
}



bar_build <- function(pl, p, DF, att){ 
	DF$tmp.data <- DF[,unlist(att[[p]]$panel.data)]
	DF$tmp.adj <- att[[p]]$graph.bar.size

	tmp.limsx <- range(c(0, DF[,unlist(att[[p]]$panel.data)], na.rm=T))
	if(any(!is.na(att[[p]]$xaxis.ticks))) tmp.limsx <- range(c(tmp.limsx, att[[p]]$xaxis.ticks))
	tmp.limsx <- tmp.limsx + c(-.001,.05) * diff(tmp.limsx)
	tmp.limsy <- -(range(DF$pGrpOrd) + c(-1,1) * .5)


	ncolors <- length(att$colors)
	  #################################
	  #################################
	  #*** insert space for median row	
	  tmp.median.limsy <- NULL

	  if(att$median.row) {
	    if(!any(DF$median)) DF <- rbind(DF, transform(DF[1,], pGrpOrd=ncolors+1, pGrp=att$m.pGrp, median=TRUE, rank='', 
								tmp.data=median(DF$tmp.data)))
	
	 
	    DF$color[DF$median] <- ncolors
	    tmp.median.limsy <- c(-.5, -1.5)
	  }

	  #################################
	  #################################


	pl  <- 
		ggplot(DF) +
		geom_rect(aes(xmin=0, ymin=-pGrpOrd-(tmp.adj/2), 
				xmax=tmp.data, ymax=-pGrpOrd+(tmp.adj/2), 
				fill=factor(color), colour='black')) +
		scale_colour_manual(values='black', guide='none') +
		scale_fill_manual(values=att$colors, guide='none') +
	     	facet_grid(pGrp~., scales="free_y", space="free")


	if(att[[p]]$median.line) pl <- pl + geom_vline(aes(xintercept = median(tmp.data)), data=DF, 
							colour=att[[p]]$median.line.col, 
							linetype = att[[p]]$median.line.typ, 
							size = att[[p]]$median.line.size)

	if(!any(is.na(att[[p]]$add.line))){
		if(length(att[[p]]$add.line.col)==1) att[[p]]$add.line.col <- rep(att[[p]]$add.line.col[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.typ)==1) att[[p]]$add.line.typ <- rep(att[[p]]$add.line.typ[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.size)==1) att[[p]]$add.line.size <- rep(att[[p]]$add.line.size[1], length(att[[p]]$add.line))

		for(j in 1:length(att[[p]]$add.line)) pl <- pl + geom_vline(xintercept = att[[p]]$add.line[j], 
										data=DF, colour=att[[p]]$add.line.col[j], 
										linetype=att[[p]]$add.line.typ[j],
										size=att[[p]]$add.line.size[j])
	  }


	pl <- plot_opts(p, pl, att)		
	pl <- graph_opts(p, pl, att)		
	pl <- axis_opts(p, pl, att, limsx=tmp.limsx, limsy=c(tmp.limsy, tmp.median.limsy), expx=FALSE)

	pl
}



bar_cl_build <- function(pl, p, DF, att){ 
	DF$tmp.data1 <- DF[,att[[p]]$panel.data[[1]]]
	DF$tmp.data2 <- DF[,att[[p]]$panel.data[[2]]]
	DF$tmp.data3 <- DF[,att[[p]]$panel.data[[3]]]

	DF$tmp.adj <- att[[p]]$graph.bar.size

	tmp.limsx <- range(c(0, DF[,unlist(att[[p]]$panel.data)], na.rm=T))
	if(any(!is.na(att[[p]]$xaxis.ticks))) tmp.limsx <- range(c(tmp.limsx, att[[p]]$xaxis.ticks))
	tmp.limsx <- tmp.limsx + c(-.001,.05) * diff(tmp.limsx)
	tmp.limsy <- -(range(DF$pGrpOrd) + c(-1,1) * .5)


	ncolors <- length(att$colors)
	  #################################
	  #################################
	  #*** insert space for median row	
	  tmp.median.limsy <- NULL

	  if(att$median.row) {
	    if(!any(DF$median)) DF <- rbind(DF, transform(DF[1,], pGrpOrd=ncolors+1, pGrp=att$m.pGrp, median=TRUE, rank='', 
								tmp.data1=median(DF$tmp.data1),
								tmp.data2=median(DF$tmp.data2),
								tmp.data3=median(DF$tmp.data3)))
	
	 
	    DF$color[DF$median] <- ncolors
	    tmp.median.limsy <- c(-.5, -1.5)
	  }

	  #################################
	  #################################

	pl  <- 
		ggplot(DF) +
		geom_rect(aes(xmin=0, ymin=-pGrpOrd-(tmp.adj/2), 
				xmax=tmp.data1, ymax=-pGrpOrd+(tmp.adj/2), 
				fill=factor(color), colour='black')) +
	     	geom_errorbarh(aes(x=tmp.data1, xmin=tmp.data2, xmax=tmp.data3, y=-pGrpOrd), 
				height=.9*att[[p]]$graph.bar.size) + 
	     	facet_grid(pGrp~., scales="free_y", space="free") +
		scale_colour_manual(values='black', guide='none') +
		scale_fill_manual(values=att$colors, guide='none')


	if(att[[p]]$median.line) pl <- pl + geom_vline(aes(xintercept = median(tmp.data)), data=DF, 
							colour=att[[p]]$median.line.col, 
							linetype = att[[p]]$median.line.typ, 
							size = att[[p]]$median.line.size)

	if(!any(is.na(att[[p]]$add.line))){
		if(length(att[[p]]$add.line.col)==1) att[[p]]$add.line.col <- rep(att[[p]]$add.line.col[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.typ)==1) att[[p]]$add.line.typ <- rep(att[[p]]$add.line.typ[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.size)==1) att[[p]]$add.line.size <- rep(att[[p]]$add.line.size[1], length(att[[p]]$add.line))

		for(j in 1:length(att[[p]]$add.line)) pl <- pl + geom_vline(xintercept = att[[p]]$add.line[j], 
										data=DF, colour=att[[p]]$add.line.col[j], 
										linetype=att[[p]]$add.line.typ[j],
										size=att[[p]]$add.line.size[j])
	  }


	pl <- plot_opts(p, pl, att)		
	pl <- graph_opts(p, pl, att)
	pl <- axis_opts(p, pl, att, limsx=tmp.limsx, limsy=c(tmp.limsy, tmp.median.limsy), expx=FALSE)

	pl
}




box_summary_build <- function(pl, p, DF, att){ 
	if(length(att[[p]]$panel.data)==5) iCols <- c(1,2,2,3,4,4,5) else iCols <- 1:7
	tmp.data <- DF[,unlist(att[[p]]$panel.data)]
	tmp.data <- tmp.data[,iCols]
	names(tmp.data) <- paste('tmp.data',1:7,sep='')
	DF <- cbind(DF, tmp.data)

	DF$tmp.adj <- att[[p]]$graph.bar.size

	tmp.limsx <- range(DF[,unlist(att[[p]]$panel.data)], na.rm=T)
	if(any(!is.na(att[[p]]$xaxis.ticks))) tmp.limsx <- range(c(tmp.limsx, att[[p]]$xaxis.ticks))
	tmp.limsx <- tmp.limsx + c(-1,1) * diff(tmp.limsx)*.05
	tmp.limsy <- -(range(DF$pGrpOrd) + c(-1,1) * .5)


	ncolors <- length(att$colors)
	  #################################
	  #################################
	  #*** insert space for median row	
	  tmp.median.limsy <- NULL


	  if(att$median.row) {
	    if(!any(DF$median)) DF <- rbind(DF, transform(DF[1,], pGrpOrd=ncolors+1, pGrp=att$m.pGrp, median=TRUE, rank='', 
								tmp.data1=median(DF$tmp.data1), 
								tmp.data2=median(DF$tmp.data2), 
								tmp.data3=median(DF$tmp.data3), 
								tmp.data4=median(DF$tmp.data4), 
								tmp.data5=median(DF$tmp.data5), 
								tmp.data6=median(DF$tmp.data6), 
								tmp.data7=median(DF$tmp.data7)))
	
	 
	    DF$color[DF$median] <- ncolors
	    tmp.median.limsy <- c(-.5, -1.5)	# c(0, -2)
	  }

	  #################################
	  #################################


	pl  <- 
		ggplot(DF) +
	     	geom_errorbarh(aes(x=tmp.data1, xmin=tmp.data1, xmax=tmp.data7, y=-pGrpOrd), 
				height=.9*att[[p]]$graph.bar.size) + 
		geom_rect(aes(xmin=tmp.data2, ymin=-pGrpOrd-(tmp.adj/4), 
					xmax=tmp.data5, ymax=-pGrpOrd+(tmp.adj/4), 
					fill=factor(color), colour='black')) +
		geom_rect(aes(xmin=tmp.data3, ymin=-pGrpOrd-(tmp.adj/2), 
					xmax=tmp.data6, ymax=-pGrpOrd+(tmp.adj/2), 
					fill=factor(color), colour='black')) +
		geom_segment(aes(x=tmp.data4, y=-pGrpOrd-(tmp.adj/2), 
					xend=tmp.data4, yend=-pGrpOrd+(tmp.adj/2)), colour='black') +
		facet_grid(pGrp~., scales="free_y", space="free") +
		scale_colour_manual(values='black', guide='none') +
		scale_fill_manual(values=att$colors, guide='none')


	if(att[[p]]$median.line) pl <- pl + geom_vline(aes(xintercept = median(tmp.data)), data=DF, 
							colour=att[[p]]$median.line.col, 
							linetype = att[[p]]$median.line.typ, 
							size = att[[p]]$median.line.size)

	if(!any(is.na(att[[p]]$add.line))){
		if(length(att[[p]]$add.line.col)==1) att[[p]]$add.line.col <- rep(att[[p]]$add.line.col[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.typ)==1) att[[p]]$add.line.typ <- rep(att[[p]]$add.line.typ[1], length(att[[p]]$add.line))
		if(length(att[[p]]$add.line.size)==1) att[[p]]$add.line.size <- rep(att[[p]]$add.line.size[1], length(att[[p]]$add.line))

		for(j in 1:length(att[[p]]$add.line)) pl <- pl + geom_vline(xintercept = att[[p]]$add.line[j], 
										data=DF, colour=att[[p]]$add.line.col[j], 
										linetype=att[[p]]$add.line.typ[j],
										size=att[[p]]$add.line.size[j])
	  }


	pl <- plot_opts(p, pl, att)		
	pl <- graph_opts(p, pl, att)		
	pl <- axis_opts(p, pl, att, limsx=tmp.limsx, limsy=c(tmp.limsy,tmp.median.limsy), expx=FALSE)

	pl
}


	