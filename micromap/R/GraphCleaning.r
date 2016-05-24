##########################
##### theme settings #####
##########################

  # creating theme functions to be called as needed from the panel building done elsewhere
  #  	each function takes a ggplot object (e.g. pl <- plots[[1]]) and adds attributes to it
  # 	as specified in the above function arguments then returns the object


### set panel background color, overall title and margins
plot_opts <- function(i, pl, a){		
	  # sets background color if specified
	if(!is.na(a[[i]]$panel.bgcolor)) pl <- pl + theme(plot.background = element_rect(colour=a[[i]]$panel.bgcolor, 
															fill=a[[i]]$panel.bgcolor))

	  # sets a header title if specified
	  # if any of the panels have a header this inserts a blank title on those that do not in order
	  # 	to keep a similar layout among panels
	if(!all(is.na(all_atts(a, 'panel.header')))){
		tmp.size <- max(as.numeric(all_atts(a, 'panel.header.size')))*10	# All headers must be same size. Default size is 10
		tmp.lineheight <- as.numeric(all_atts(a, 'panel.header.lineheight'))
		tmp.lineheight <- tmp.lineheight[which.max(abs(tmp.lineheight-1))] 	# All line heights must be equal so we use the height 
											# most drastically different from 1
		

		  # If all headers aren't NA then we must change the titles to be blank in order to leave space
		  # 	for them at the top. If any header is multiple lines (ie contains '\n') then we must add in 
		  # 	the correct number of character returns to the other headers in order to make the plot uniform
		tmp.headers <- lapply(all_atts(a, 'panel.header'), function(t) if(is.na(t)|t=='') t=' ' else t)
		tmp.title <- tmp.headers[[i]]
		ns <- max(unlist(lapply(tmp.headers, function(x) length(strsplit(x, '\n')[[1]])))) - length(strsplit(tmp.title,'\n')[[1]])  
		if(ns>0) tmp.title <- paste(tmp.title, rep(' \n ',ns), sep='')
		

		pl <- pl + ggtitle(tmp.title) +
				theme(plot.title=element_text(family=a[[i]]$panel.header.font, face=a[[i]]$panel.header.face, 
							colour=a[[i]]$panel.header.color, size=tmp.size, lineheight=tmp.lineheight))
	} 
	

	  # sets panel margins and removes ggplot's default side strips
	pl <- pl + theme(strip.background = element_blank(), 
				strip.text.x = element_blank(), 
				strip.text.y = element_blank(), 
				plot.margin = unit(a[[i]]$panel.margins, "lines"))

	pl
}




### set graph background color and whether grid lines show up
graph_opts <- function(i, pl, a){		
	bgcolor <- ifelse(!is.na(a[[i]]$panel.bgcolor), a[[i]]$panel.bgcolor, 'white')
	bgcolor <- ifelse(!is.na(a[[i]]$graph.bgcolor), a[[i]]$graph.bgcolor, bgcolor)

	  # sets background color of graphs 
	  # note: what were referring to as "graphs" are what ggplot refers to as "panels" (ie. "panel.background")
	pl <- pl + theme(panel.background = element_rect(colour=bgcolor, fill=bgcolor))

	
	  # draws grid lines in the specified color if desired -- defaults to darkgray in attribute list
	if(a[[i]]$graph.grid.major){ 
			pl <- pl + theme(panel.grid.major = element_line(colour=a[[i]]$graph.grid.color))
		} else {
			pl <- pl + theme(panel.grid.major = element_blank())
	  }
	if(a[[i]]$graph.grid.minor){ 
			pl <- pl + theme(panel.grid.minor = element_line(colour=a[[i]]$graph.grid.color))
		} else {
			pl <- pl + theme(panel.grid.minor = element_blank())
	  }
	
	pl
}



### sets graph boundaries, ticks, labels, borders
axis_opts <- function(i, pl, a, limsx=NA, limsy=NA, border=TRUE, expx=FALSE){
	# i=p; a=att; limsx=tmp.limsx; limsy=c(tmp.limsy,tmp.median.limsy); border=FALSE; expx=FALSE

	# many features are "hidden" by simply coloring the same color as the background so
	#   if panel background is NA we assume "white" will effectively do the hiding
    bgcolor <- ifelse(!is.na(a[[i]]$panel.bgcolor), a[[i]]$panel.bgcolor, 'white')

	# specify label size as maximum of all requested label sizes
    label.size <- as.numeric(max(all_atts(a, 'xaxis.labels.size')))*10


      # limsy will sometimes be in the form (c(lower bound, upper bound, lower bound for the median , upper bound for the median))
	#   if thats the case, we split it into the two seperate limits here
    median.limsy <- NULL
    if(length(limsy)==4) {median.limsy <- limsy[3:4]; limsy <- limsy[1:2]}


    ##############    
    ### X axis ### 
    ############## 
	  # with ggplot2, most axis specifications need to be made through the "scale_x_continuous()" function. there
	  #   are 5 arguements: title, breaks (tick locations), labels (tick labels), limits (data limits), and expand (the 
	  #   extent to which the axis is expanded beyond the limits of the data). these specifications must be made 
	  # 	all at once so we build this statement as a string and then execute it through an "eval" statement at the end

	  # the following boolean variables state whether to include these specifications in the scale_x_continuous statement
	  # 	we start by assuming none except title will be needed 
	x.breaks <- x.labels <- x.limits <- x.expand <- FALSE


	 ### axis title ###
	xstr.title <- "''"

	  # If all axis titles aren't NA then we must change the other titles to be blank in order to leave space
	  # 	for them at the bottom. If any title is multiple lines (ie contains '\n') then we must add in 
	  # 	the correct number of character returns to the other titles in order to make the plot uniform
	if(!all(is.na(all_atts(a, 'xaxis.title')))){
		tmp.size <- max(as.numeric(all_atts(a, 'xaxis.title.size')))*8		# All titles must be same size. Default size is 8
		tmp.lineheight <- as.numeric(all_atts(a, 'xaxis.title.lineheight'))
		tmp.lineheight <- tmp.lineheight[which.max(abs(tmp.lineheight-1))] 	# All line heights must be equal so we use the height 
											# most drastically different from 1


		tmp.titles <- lapply(all_atts(a, 'xaxis.title'), function(t) if(is.na(t)|t=='') t=' ' else t)
		tmp.title <- tmp.titles[[i]]
		ns <- max(unlist(lapply(tmp.titles, function(x) length(strsplit(x, '\n')[[1]])))) - length(strsplit(tmp.title,'\n')[[1]])  
		if(ns>0) tmp.title <- paste(tmp.title, rep(' \n ',ns), sep='')
		
		xstr.title <- paste("'",tmp.title,"'",sep='')
		pl <- pl + theme(axis.title.x = element_text(family=a[[i]]$xaxis.title.font, face=a[[i]]$xaxis.title.face, 
							colour=a[[i]]$xaxis.title.color, size=tmp.size, lineheight=tmp.lineheight))
	} 
		

	 ### axis limits and expansion ###
	if (!any(is.na(limsx))) x.limits <- TRUE

		# if there is a border to be added, we must manually deal with expansion
        if(!expx){				
		x.expand <- TRUE
		xstr.expand <- as.character(", expand=c(0,0)")
	  }

	xstr.limits <- as.character(paste('c(',min(limsx), ',', max(limsx),')'))
	xstr.limits <- paste(", limits=", xstr.limits)



	 ### panel footers (not completed) ###
	  # "panel footers" are really just augmented x axis titles
	  # if all axis titles are blank then we hide axis titles on the whole plot	
	if(all(is.na(all_atts(a, 'xaxis.title')))) pl <- pl + theme(axis.title.x = element_blank()) 	


	 ### axis lines ###
	  # note: axis lines are always there, if the user doesn't want to 
	  #       see them they are colored to match the background
	if (!a[[i]]$xaxis.line.display & !a[[i]]$yaxis.line.display) {
		pl <- pl + theme(axis.line = element_line(colour=bgcolor)) 
		
	     # else lines will be plotted
	   } else {
		pl <- pl + theme(axis.line = element_line(colour='black')) 
	 } 


	 ### axis ticks ###
	  # for now we assume ticks are never wanted as they make things pretty cluttered looking
	pl <- pl + theme(axis.ticks = element_blank())


	 ### axis text ###
	  # trys to hide axis text on the whole plot
	if(!any(all_attsb(a, 'xaxis.text.display'))) {									
		pl <- pl + theme(axis.text.x = element_blank())

	    # otherwise trys to "hide" axis text on this panel
	  } else if (!a[[i]]$xaxis.text.display) {								
		pl <- pl + theme(axis.text.x = element_text(colour=bgcolor, size=label.size))	

	    # axis text will show and we'll add specific labels if requested
	  } else if (!is.na(unlist(a[[i]]$xaxis.labels)[1]) & 
				!is.na(unlist(a[[i]]$xaxis.ticks)[1])) { 	
	
	 	tmpTheme <- "theme(axis.text.x = element_text(size=label.size"
		if(!is.null(a[[i]]$xaxis.labels.angle)) tmpTheme <- paste(tmpTheme, ", angle = ", a[[i]]$xaxis.labels.angle)
		if(!is.null(a[[i]]$xaxis.labels.hjust)) tmpTheme <- paste(tmpTheme, ", hjust =", a[[i]]$xaxis.labels.hjust)
		if(!is.null(a[[i]]$xaxis.labels.vjust)) tmpTheme <- paste(tmpTheme, ", vjust =", a[[i]]$xaxis.labels.vjust)
		tmpTheme <- paste(tmpTheme, "))")

		pl <- pl + eval(parse(text=tmpTheme	))
		
		x.breaks <- x.labels <- TRUE		
		xstr.breaks <- paste(', breaks=c(', make.string(a[[i]]$xaxis.ticks),')',sep='')
		xstr.labels <- paste(', labels=c(', make.string(a[[i]]$xaxis.labels),')',sep='')

	    # warning if user only specified text but not location or vice versa
	  } else if (!is.na(unlist(a[[i]]$xaxis.labels)[1]) | 
			!is.na(unlist(a[[i]]$xaxis.ticks)[1])) { 
		print('Warning: both axis labels AND tick location must be specified')

	   # otherwise text shows up as ggplot defaults
	}
 
 
	  # put it all together and execute the eval call
	xstr <- paste("scale_x_continuous(", xstr.title)
	if (x.expand) xstr <- paste(xstr, xstr.expand)
	if (x.breaks) xstr <- paste(xstr, xstr.breaks)
	if (x.labels) xstr <- paste(xstr, xstr.labels)
	if (x.limits) xstr <- paste(xstr, xstr.limits)
	xstr <- paste(xstr, ")")

	pl <- pl + eval(parse(text=xstr))
 
 
 	
    ##############
    ### Y axis ###   
    ##############
	  # with ggplot2, most axis specifications need to be made through the "scale_y_continuous()" function. there
	  #   are 5 arguements: title, breaks (tick locations), labels (tick labels), limits (data limits), and expand (the 
	  #   extent to which the axis is expanded beyond the limits of the data). these specifications must be made 
	  #   all at once so we build this statement as a string and then execute it through an "eval" statement at the end

	 ### axis title ###
	ystr.title  <- ifelse(!is.na(a[[i]]$yaxis.title), a[[i]]$yaxis.title, "''")
	

	 ### axis text ###
	if(a[[i]]$yaxis.ticks.display | a[[i]]$yaxis.text.display) ystr.breaks <- "" else ystr.breaks <- ", breaks=NULL" 


	 ### axis limits and expansion ###
	ystr.expand <- ", expand=c(0,0)"
	limsy <- limsy + c(-1,1) * diff(limsy)*a$plot.pGrp.spacing

	ystr.limits <- as.character(paste('c(',min(limsy), ',', max(limsy),')'))
	ystr.limits <- paste(", limits=", ystr.limits)

	pl <- pl + theme(panel.margin = unit(0, "lines"))

	
	# if (any(is.na(limsy)) | a$median.row) y.limits <- FALSE else y.limits <- TRUE


	  # put it all together and execute the eval call
	ystr <- paste("scale_y_continuous(", ystr.title, ystr.expand, ystr.breaks)
	# if (y.limits) ystr <- paste(ystr, ystr.limits)
	ystr <- paste(ystr, ")")

	pl <- pl + eval(parse(text=ystr))


    ##############
    ### border ###   
    ##############
 
  	borderx <- range(limsx) + c(1,-1) * diff(range(limsx))*.001
	bordery <- range(limsy) + c(0, -1) * diff(range(limsy))*.001
	if(!is.null(median.limsy)) median.limsy <- range(median.limsy) - c(0, diff(range(median.limsy))*.001)

	tmp.border <- data.frame(pGrp=rep(1:max(pl$data$pGrp),each=2), ymin=bordery[1], ymax=bordery[2], 
									xmin=borderx[1], xmax=borderx[2])
	if(a$median.row) tmp.border <- rbind(subset(tmp.border, !pGrp==a$m.pGrp), data.frame(pGrp=a$m.pGrp, ymin=median.limsy[1], ymax=median.limsy[2],
											xmin=borderx[1], xmax=borderx[2]))

	if(border) border.color <- a[[i]]$graph.border.color else border.color <- NA
	pl <- pl + geom_rect(aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax), data=tmp.border, 
					colour= border.color, fill=NA)
	pl <- pl + theme(axis.line = element_blank())

	pl

}





		
assimilatePlot <- function(pl, i, a, limsx=NA, limsy=NA){
	pl <- plot_opts(i,pl,a)
	pl <- graph_opts(i,pl,a)

	if(is.na(limsx)){ 
		limsx <- range(pl$data[,unlist(a[[i]]$panel.data)])
		limsx <- limsx+c(-1,1)*diff(limsx)*.05
	  }

	if(is.na(limsy)){
		# labs <- names(pl$options$labels)
		# labs <- labs[sapply(1:length(labs), function(j) pmatch("y",labs[j], nomatch=0)==1)]

		# limsy <- -range(pl$data[,sub("-","",unlist(pl$options$labels[labs]))])
		# limsy <- limsy + c(1,-1)*diff(limsy)*.05
		limsy <- -c(.5, max(a$grouping)+.5)
	  }

	if(a$median.row){	
		pl <- pl + scale_colour_manual(values=c(a$colors,gray(.5)), guide='none') 
		
		median.limsy <-  c(-.5, -1.5)
		limsy <- c(limsy, median.limsy)
	}

	pl <- axis_opts(i,pl,a, limsx=limsx, limsy=limsy, border=TRUE)
	pl
}


