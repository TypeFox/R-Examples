RankMaps <- function(pl, p, mapDF, att){
	# att=a
	bgcolor 	<- ifelse(!is.na(att[[p]]$panel.bgcolor), att[[p]]$panel.bgcolor, 'white')
	outer.hull 	<- att[[p]]$outer.hull		# whether we need to construct an outline of the map poygons
	ncolors 	<- length(att$colors)
	nGroups 	<- range(mapDF$pGrp[!is.na(mapDF$pGrp)])	

	mapDF$fill.color <- mapDF$pGrpOrd
	mapDF$IDpoly 	 <- paste(mapDF$ID,mapDF$poly,sep='.')

	map.all <- att[[p]]$map.all

	#*** store all the polygon data then subset only those with data
	mapDF.all <- mapDF				
        mapDF <- subset(mapDF.all, !is.na(rank))	

	#*** here we build the "outer hull" (polygon outline) if neccesary
	if(outer.hull==TRUE){
		ii=0
		lHull <- NULL
		for(rr in unique(mapDF.all$region)){
		  ii <- ii + 1
	
		  tmpPlys <- vector("list", length(unique(subset(mapDF.all, region==rr)$poly)))
		  jj=1
		  for(pp in unique(subset(mapDF.all, region==rr)$poly)){
			tmpPlys[[jj]] <- Polygon(subset(mapDF.all, region==rr & poly==pp)[,c('coordsx','coordsy')])
		
			jj <- jj + 1	
		  }
				
		  lHull <- append(lHull, Polygons(tmpPlys, ID=rr))
		}
		spHull <- SpatialPolygons(lHull)
		
		unionHull <- unionSpatialPolygons(spHull, rep(1, length(spHull)))
		
		    lHull <- vector("list", length(unionHull@polygons))
		    for (i in 1:length(lHull)) {
		        tmp <- unionHull@polygons[[i]]
		        lHull[[i]] <- lapply(1:length(tmp@Polygons), function(j) cbind(i, 
		            j, tmp@Polygons[[j]]@labpt[1], tmp@Polygons[[j]]@labpt[2], 
		            tmp@Polygons[[j]]@coords, tmp@Polygons[[j]]@hole, 
		            tmp@Polygons[[j]]@area))
		    }

		    dHull <- NULL
		    for (i in 1:length(lHull)) {
		        for (j in 1:length(lHull[[i]])) {
		            dHull <- rbind(dHull, lHull[[i]][[j]])
		        }
		    }
		
		    dHull <- data.frame(1, dHull)
		    names(dHull) <- c("ID", "region", "poly", "lab.x", "lab.y", "coordsx", "coordsy", "hole", "area")
		    dHull <- transform(dHull, poly = (region - 1) * max(dHull$poly) + poly)
		
			dHull$fill.color <- dHull$pGrpOrd
			dHull$IDpoly <- paste(dHull$ID,dHull$poly,sep='.')
	}


	#*** If there is an odd number of polygons and a median row in our map 
	#***   we must deal with the median polygon specially
	if(att$median.row & any(mapDF$pGrp==att$m.pGrp)){
		#*** find the median polygon set the fill color for it to be the user specified median color
		#*** create copies of that polygon to be displayed in the preceding 
		#*** and subsequent perceptual groups  
		mapDF.median 		<- subset(mapDF, pGrp==att$m.pGrp)
		mapDF.median$fill.color <- ncolors+1	
		mapDF.median		<- rbind(transform(mapDF.median, pGrp=(att$m.pGrp-1), IDpoly='median1'),
						transform(mapDF.median, pGrp=(att$m.pGrp+1), IDpoly='median2'))
		
		mapDF <- rbind(subset(mapDF, !pGrp==att$m.pGrp), mapDF.median)
	}

	#*** Creates the map panel object 	
	pl <- ggplot(mapDF, aes(x=coordsx, y=coordsy, group=IDpoly)) 


	#*** If we are displaying polygons without associated data we 
	#***   do so first so they appear in the background
	if(map.all==TRUE){
		if(nrow(subset(mapDF.all, is.na(pGrp)))>0) for(g in nGroups[1]:nGroups[2]) pl <- pl + 
			geom_polygon(fill=att[[p]]$nodata.fill, 
				colour=att[[p]]$nodata.border.color, 
				size=att[[p]]$nodata.border.size/2, 
				data=transform(subset(mapDF.all, is.na(pGrp)), pGrp=g))
	}


	#*** Draw polygons from all previous perceptual groups by looping through
	#***   each perceptual group and creating datasets which contain all 
	#***   previous perceptual groups or future groups depending on the 
	#***   user's shading preference
	if(att[[p]]$fill.regions=="aggregate"){
		lWith.data <- 1:(nGroups[2]-1)
		lWith.data <- lWith.data[!lWith.data==att$m.pGrp]
		lWithout.data <- max(nGroups[1],2):nGroups[2]
		lWithout.data <- lWithout.data[!lWithout.data==att$m.pGrp]

		for(g in lWith.data) pl <- pl + 
			geom_polygon(fill=att[[p]]$withdata.fill, 
				colour=att[[p]]$withdata.border.color, 
				size=att[[p]]$withdata.border.size/2, 
				data=transform(mapDF[mapDF$pGrp>g,], pGrp=g))						

		for(g in lWithout.data) pl <- pl + 
			geom_polygon(fill=att[[p]]$inactive.fill, 
				colour=att[[p]]$inactive.border.color,  
				size=att[[p]]$inactive.border.size/2, 
				data=transform(mapDF[mapDF$pGrp<g,], pGrp=g)) 

	 }
 	
	if(att[[p]]$fill.regions=="two ended"){
		lGroups <- nGroups[1]:floor(att$m.pGrp)
		lGroups <- lGroups[!lGroups==att$m.pGrp]
		for(g in lGroups) pl <- pl + 
			geom_polygon(fill=att[[p]]$withdata.fill, 
				colour=att[[p]]$withdata.border.color,  
				size=att[[p]]$withdata.border.size/2, 
				data=transform(mapDF[mapDF$pGrp>g,], pGrp=g)) 	

		lGroups <- max(nGroups[1],2):floor(att$m.pGrp)
		lGroups <- lGroups[!lGroups==att$m.pGrp]
		for(g in lGroups) pl <- pl + 
			geom_polygon(fill=att[[p]]$inactive.fill, 
				colour=att[[p]]$inactive.border.color,  
				size=att[[p]]$inactive.border.size/2, 
				data=transform(mapDF[mapDF$pGrp<g,], pGrp=g)) 	

		lGroups <- ceiling(att$m.pGrp):nGroups[2]
		lGroups <- lGroups[!lGroups==att$m.pGrp]
		for(g in lGroups) pl <- pl + 
			geom_polygon(fill=att[[p]]$withdata.fill, 
				colour=att[[p]]$withdata.border.color,  
				size=att[[p]]$withdata.border.size/2, 
				data=transform(mapDF[mapDF$pGrp<g,], pGrp=g)) 

		lGroups <- ceiling(att$m.pGrp):(nGroups[2]-1)
		lGroups <- lGroups[!lGroups==att$m.pGrp]
		for(g in lGroups) pl <- pl + 
			geom_polygon(fill=att[[p]]$inactive.fill, 
				colour=att[[p]]$inactive.border.color,  
				size=att[[p]]$inactive.border.size/2, 
				data=transform(mapDF[mapDF$pGrp>g,], pGrp=g)) 

	 
	}

	if(att[[p]]$fill.regions=="with data"){
		lWith.data <- nGroups[1]:nGroups[2]
		lWith.data <- lWith.data[!lWith.data==att$m.pGrp]

		for(g in lWith.data) pl <- pl + 
			geom_polygon(fill=att[[p]]$withdata.fill, 
				colour=att[[p]]$withdata.border.color,  
				size=att[[p]]$withdata.border.size/2, 
				data=transform(mapDF[!mapDF$pGrp==g,], pGrp=g)) 

	 }



	#*** draw polygons of current perceptual group																			
	pl <- pl + 
		geom_polygon(aes(fill=factor(fill.color)), 
			colour=att[[p]]$active.border.color, 
			size=att[[p]]$active.border.size/2, 
			data=subset(mapDF, hole==0)) + 				
		facet_grid(pGrp~., scales="free_y", space="free") +
		scale_fill_manual(values=c(att$colors), guide='none') +
		coord_equal() 


	  #################################
	  #################################
	  #*** insert white space for median row	
	  mapDF.median <- data.frame(pGrpOrd=1, pGrp=att$m.pGrp, 
					rank=(max(mapDF$rank)+1)/2, ID='median', 
					coordsx=range(mapDF$coordsx)[c(1,1,2,2,1)], 
					coordsy=c(.5, -.5)[c(1,2,2,1,1)],
#					textx=median(range(mapDF$coordsx)), texty=median(range(mapDF$coordsy)),	
					tmp.label=att$median.text.label,
					textx=median(range(mapDF$coordsx)), texty=0, tmp.label='Median',
					region=1, poly=1, plug=0, hole=0, IDpoly='median')
	
	  if(att$median.row) pl <- pl + 
			geom_polygon(fill='white', colour='white', data=mapDF.median) + 
				geom_text(aes(x=textx, y=texty, label=tmp.label, hjust=.5, vjust=.4),
						colour=att$median.text.color, size=5*att$median.text.size, data=mapDF.median) +
				facet_grid(pGrp~., scales="free_y", space="free")

	  #################################
	  #################################

	  #*** Throw an outer hull over the outside

	if(outer.hull==TRUE) for(g in nGroups[1]:nGroups[2]) pl <- pl + geom_polygon(fill=NA,	 
													colour=att[[p]]$outer.hull.color, 
													size=att[[p]]$outer.hull.size, 
													data=data.frame(dHull, pGrp=g))


	#*** Creates backgroupnd points to correctly align perceptual groups 		
	bdrCoordsy <- att[[p]]$bdrCoordsy
	bdrCoordsy <- bdrCoordsy + c(1,-1) * diff(bdrCoordsy)*att$plot.pGrp.spacing
	bdrCoordsy <- bdrCoordsy + c(0, -1) * diff(range(bdrCoordsy))*.001

	bdrCoordsx <- att[[p]]$bdrCoordsx

	bdrCoords <- data.frame(coordsx=bdrCoordsx, coordsy=bdrCoordsy, pGrp=NA, IDpoly=NA)
	lGroups <- nGroups[1]:nGroups[2]
	lGroups <- lGroups[!lGroups==att$m.pGrp]
	for(g in lGroups) pl <- pl + 
		geom_point(fill=bgcolor, size=.001,
				colour=bgcolor, 			
				data=transform(bdrCoords, pGrp=g))


	  #*** "house keeping" functions to align all panels and apply the user's graphical specifications\	
	xstr.title <-  "''"


	  # If all axis titles aren't NA then we must change the other titles to be blank in order to leave space
	  # 	for them at the bottom. If any title is multiple lines (ie contains '\n') then we must add in 
	  # 	the correct number of character returns to the other titles in order to make the plot uniform
	if(!all(is.na(all_atts(att, 'xaxis.title')))){
		tmp.titles <- lapply(all_atts(att, 'xaxis.title'), function(t) if(is.na(t)|t=='') t=' ' else t)
		tmp.title <- tmp.titles[[p]]
		ns <- max(unlist(lapply(tmp.titles, function(x) length(strsplit(x, '\n')[[1]])))) - length(strsplit(tmp.title,'\n')[[1]])  
		if(ns>0) tmp.title <- paste(tmp.title, rep(' \n ',ns), sep='')
		
		xstr.title <- paste("'",tmp.title,"'",sep='')
		pl <- pl + theme(axis.title.x = element_text(size=8*att[[p]]$xaxis.title.size, colour=bgcolor)) 
	  } else {
		pl <- pl + theme(axis.title.x = element_blank()) 
	}
	

	if(any(c(all_attsb(att, 'xaxis.ticks.display'), all_attsb(att, 'xaxis.text.display')))){
		pl <- pl + theme(axis.text.x = element_text(colour=bgcolor)) +
				theme(axis.ticks = element_line(colour=bgcolor))
	 } else { 
		pl <- pl + theme(axis.text.x = element_blank()) +
				theme(axis.ticks = element_blank()) 
	 }



	ystr <- paste("scale_y_continuous('', breaks=NULL, expand=c(0,0))")

	pl <- pl + eval(parse(text=ystr))

	pl  <- plot_opts(p, pl, att)		
	pl  <- graph_opts(p, pl, att)	

	pl <- pl + theme(panel.margin = unit(0, "lines"))


	pl 

}


CatMaps <- function(pl, p, mapDF, att){
	bgcolor <- ifelse(!is.na(att[[p]]$panel.bgcolor), att[[p]]$panel.bgcolor, 'white')

	mapDF.all <- mapDF
	mapDF <- subset(mapDF, !is.na(pGrp))

	pl <-  
		ggplot(mapDF, aes(x=coordsx, y=coordsy, group=poly)) +  
		geom_polygon(fill='white', colour='black', data=transform(mapDF.all, pGrp=NULL)) + 	
		geom_polygon(fill=att$map.color, colour='black', data=subset(mapDF, hole==0)) + 				
		geom_polygon(fill='white', colour='black', data=subset(mapDF, hole==1)) +
		geom_polygon(fill='transparent', colour='black', data=subset(transform(mapDF, pGrp=NULL), hole==1)) +
		facet_grid(pGrp~., scales="free_y", space="free") +
		coord_equal() 



	#*** Creates backgroupnd points to correctly align perceptual groups 		
	bdrCoordsy <- att[[p]]$bdrCoordsy
	bdrCoordsy <- bdrCoordsy + c(1,-1) * diff(bdrCoordsy)*att$plot.pGrp.spacing
	bdrCoordsy <- bdrCoordsy + c(0, -1) * diff(range(bdrCoordsy))*.001

	bdrCoordsx <- att[[p]]$bdrCoordsx

	bdrCoords <- data.frame(coordsx=bdrCoordsx, coordsy=bdrCoordsy, pGrp=NA, poly=NA)
	pl <- pl + geom_point(fill=bgcolor, size=.001,
				colour=bgcolor, 			
				data=transform(bdrCoords, pGrp=NULL))


	  #*** "house keeping" functions to align all panels and apply the user's graphical specifications\	
	xstr.title <-  "''"


	  # If all axis titles aren't NA then we must change the other titles to be blank in order to leave space
	  # 	for them at the bottom. If any title is multiple lines (ie contains '\n') then we must add in 
	  # 	the correct number of character returns to the other titles in order to make the plot uniform
	if(!all(is.na(all_atts(att, 'xaxis.title')))){
		tmp.titles <- lapply(all_atts(att, 'xaxis.title'), function(t) if(is.na(t)|t=='') t=' ' else t)
		tmp.title <- tmp.titles[[p]]
		ns <- max(unlist(lapply(tmp.titles, function(x) length(strsplit(x, '\n')[[1]])))) - length(strsplit(tmp.title,'\n')[[1]])  
		if(ns>0) tmp.title <- paste(tmp.title, rep(' \n ',ns), sep='')
		
		xstr.title <- paste("'",tmp.title,"'",sep='')
		pl <- pl + theme(axis.title.x = element_text(size=8*att[[p]]$xaxis.title.size, colour=bgcolor)) 
	  } else {
		pl <- pl + theme(axis.title.x = element_blank()) 
	}
	

	if(any(c(all_attsb(att, 'xaxis.ticks.display'), all_attsb(att, 'xaxis.text.display')))){
		pl <- pl + theme(axis.text.x = element_text(colour=bgcolor)) +
				theme(axis.ticks = element_line(colour=bgcolor))
	 } else { 
		pl <- pl + theme(axis.text.x = element_blank()) +
				theme(axis.ticks = element_blank()) 
	 }


#	ystr <- paste("scale_y_continuous('', breaks=NULL, expand=c(",pGrp.spacing,",",pGrp.spacing,"))")
	ystr <- paste("scale_y_continuous('', breaks=NULL, expand=c(0,0))")
	pl <- pl + eval(parse(text=ystr))

	pl  <- plot_opts(p, pl, att)		
	pl  <- graph_opts(p, pl, att)	
	pl <- pl + theme(panel.margin = unit(0, "lines")) 

	pl

}


