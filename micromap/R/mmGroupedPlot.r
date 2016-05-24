# ***************************************** #
# *** variable names are mostly preceded by a character indicating the data type they are designed to hold
# *** 't*' for a text variable, 'i*' for a variable containing integers (although this is used for numeric entries in general)
# *** 'b*' for boolean, 'l*' for lists of text or numeric strings
# ***************************************** #

mmgroupedplot <- function(stat.data, map.data, 		# Required -- statistical data; map data
  panel.types, panel.data, 				# Required -- panel types (e.g. 'map', 'labels', 'dot.cl', etc.);
								# 	which columns in dStats to get the data from
  map.link=NULL,						# Specify in a vector (e.g. "c('Subpopulation', 'Poly_Name')") the columns from
								# 	dStats and dMap (respectively) to link the correct map to the correct data line
  nPanels=length(panel.types),				# A count of the number of panels (not meant to actually be specified by a user
								# 	but is used referenced in this function arguments list
  grp.by, cat,						# Required	-- specifies which column in dStats by which to rank and order the output;
								# 	specifies the perceptual groups (e.g. "c(5,5,5,5,2)" to have 4 gorups of 5 followed 
								# 	by a single group of 2			
  colors=brewer.pal(10, "Spectral"),	
  map.color='lightyellow',
  map.all=FALSE,	

  print.file='no', print.res=NA,

  panel.att=vector("list", nPanels),

  ### Options for entire linked micromap plot ###
  plot.header=NA,
  plot.header.size=NA,
  plot.header.color=NA,
  plot.footer=NA,
  plot.footer.size=NA,
  plot.footer.color=NA,
	
  plot.width=7,
  plot.height=7,

  map.spacing=1,
  plot.grp.spacing=1,
  plot.panel.spacing=1,
  plot.panel.margins=c(0,0,1,0)
){					### end function arguments list

  # rename function inputs (for ease in coding)
dStats <- stat.data
dMap <- map.data


tPlot.header=plot.header
iPlot.header.size=plot.header.size
tPlot.header.color=plot.header.color
tPlot.footer=plot.footer
iPlot.footer.size=plot.footer.size
tPlot.footer.color=plot.footer.color


pps = plot.panel.spacing*.05


##############################
## create attribute object ###
##############################

  # Set "plot wide" options first
plot.atts <- list(
  			plot.grp.spacing=plot.grp.spacing*.05,
			plot.pGrp.spacing=plot.grp.spacing*.05,
  			plot.panel.margins=plot.panel.margins,
  			grp.by=grp.by, 
  			colors=colors,
			map.color=map.color,
			map.spacing=map.spacing*.025,
			plot.width=plot.width,
			plot.height=plot.height)

  # grab default attribute lists for each panel
a <- vector("list", nPanels)
for(p in 1:nPanels) {	
    att.name <- paste(panel.types[p],'_att',sep='')				# <panel type>.att is the name of the function that creates the 
    att.function <- paste(panel.types[p],'_att()',sep='')			# 	default attribute list. we check that it exists then run it
    if(exists(att.name)) a[[p]] <- eval(parse(text=att.function)) else a[[p]] <- eval(parse(text=standard_att))
}
a <- append(a, plot.atts)



  # panel data must be specified for every panel type so 
  #	we acount for it here first
for(j in 1:length(panel.types)) {
	a[[j]] <- append(a[[j]], list(panel.data=unlist(panel.data[[j]])))
 }


# *** loop through attributes given in input, match the names with those in "a"
# ***  and change the values within
  for(j in 1:length(panel.att)){		
    k <- unlist(panel.att[[j]][1], use.names=FALSE)

    for(s in names(panel.att[[j]])[-1]){	# s='header'
   	w <- match(s, right(names(a[[k]]), nchar(s)))
	if(is.na(w)) w <- match(paste('panel.',s,sep=''), names(a[[k]]))	# the "panel." part of the attribute name is somewhat 
												#	implies so a user may have left it off 

	  # replace the attribute or warn that no attribute by that name was found										
	if(!is.na(w)) a[[k]][[w]] <- unlist(panel.att[[j]][s], use.names=FALSE) else print(paste('attribute:',s,'is not recognized for panel',j))

    if(!is.null(k)){
    	w <- match('left.margin', names(a[[k]]))
    	if(!is.na(w) & !is.na(a[[k]][[w]])) 	a[[k]]$panel.margins[4] <- a[[k]]$panel.margins[4] + a[[k]][[w]]
    	w <- match('right.margin', names(a[[k]]))
    	if(!is.na(w) & !is.na(a[[k]][[w]])) 	a[[k]]$panel.margins[2] <- a[[k]]$panel.margins[2] + a[[k]][[w]]
      }

    }
  }


##########################
## data reorganization ###
##########################

# *** move dStats into the DF variable, add an overall rank column, a perceptual group columnbased on ranks, and a within 
# *** and a within perceptual group ranking

DF <- create_DF_cat(dStats, grp.by, cat)

# *** these are all meaningless in the category context but are referenced by other functions
a$median.row <- FALSE			
a$two.ended.maps <- FALSE
a$m.pGrp <- (max(DF$pGrp)+1)/2
a$m.rank <- NA


# Many panel plotting functions add extra columns to the DF table. This is
#   created here as a reference to default back to after each panel is conctructed
DF.hold <- DF

if(any(panel.types=='map')){
  w <- match(dMap[,map.link[2]], unique(DF[,map.link[1]]))
  if(!map.all) mapDF <- dMap[!is.na(w),] else mapDF <- dMap


  # make sure there is a hole and plug column. If not, just insert dummy columns 
  if(!'hole'%in%names(mapDF)) mapDF$hole <- 0
  if(!'plug'%in%names(mapDF)) mapDF$plug <- 0


  tmpDF <- unique(DF[,c('pGrp', 'pGrpOrd', map.link[1])])
  w <- match(mapDF[,map.link[2]], tmpDF[,map.link[1]])

  mapDF <- cbind(pGrp=DF[w,'pGrp'], mapDF)
  mapDF <- cbind(pGrpOrd=DF$pGrpOrd[w], mapDF)

  mapPanelWidth <- as.numeric(a[[which(panel.types=='map')]]$panel.width)
  totalUnitWidth = pps * (length(all_atts(a,'panel.width'))+1) + sum(as.numeric(all_atts(a,'panel.width')))


  mapPanelHeight <- max(DF$pGrpOrd[!is.na(DF$pGrpOrd)])+.5
  totalUnitHeight = max(DF$pGrp[!is.na(DF$pGrp)]) * mapPanelHeight 

  ns <- max(unlist(lapply(all_atts(a, 'panel.header'), function(x) length(strsplit(x, '\n')[[1]])))) 
  adj.plot.height <- plot.height - 
		ns*.2 + 							# header lines
		any(!is.na(all_atts(a, 'xaxis.title')))*.1 +
		.2 + .5							# top & bottom margin

  plot.h2w.ratio <- (mapPanelHeight/totalUnitHeight * adj.plot.height) / (mapPanelWidth/totalUnitWidth * plot.width)


    # change coordinates to align properly with other panels
  nYrows <- diff(range(mapDF$pGrpOrd[!is.na(mapDF$pGrpOrd)]) + c(-1,1) * .5)  
  nXcols <- nYrows/plot.h2w.ratio		

  limitingAxis <- ifelse(diff(range(mapDF$coordsy)) / diff(range(mapDF$coordsx)) > 
					 plot.h2w.ratio, 'y','x')


  if(limitingAxis=='y'){
	  redFact <- diff(range(mapDF$coordsy))/nYrows 
	  mapDF$coordsy <- (mapDF$coordsy - median(range(mapDF$coordsy)))/diff(range(mapDF$coordsy))*nYrows
	  mapDF$coordsx <-   (mapDF$coordsx - median(range(mapDF$coordsx)))/redFact 
  }
  if(limitingAxis=='x'){
	  redFact <- diff(range(mapDF$coordsx))/nXcols 
	  mapDF$coordsx <- (mapDF$coordsx - median(range(mapDF$coordsx)))/diff(range(mapDF$coordsx))*nXcols
	  mapDF$coordsy <-   (mapDF$coordsy - median(range(mapDF$coordsy)))/redFact 
  }
  
  a[[which(panel.types=='map')]]$bdrCoordsy <- c(-1,1) * nYrows/2
  a[[which(panel.types=='map')]]$bdrCoordsx <- c(-1,1) * nXcols/2
  	 
  
  rm(tmpDF, dMap)
}



# *** set up a list to store plot objects to be created later
# *** note: each panel in the plot is a "plot object" 
plots <- vector("list", nPanels)



###############################
##### create plot objects #####
###############################
# *** each section builds a plot object of the specified type
for(p in 1:nPanels){		

  if (panel.types[p]=='map'){

	plots[[p]]  <- CatMaps(plots[[p]], p, mapDF, a) 
  
  } else if(exists(as.character(paste(panel.types[p],'_build',sep='')))) {					# all graph types should have a function by the 
												# same name. First we check to see if such a function does
												# in fact exist, if so we use "eval(parse(..." to call it
	plots[[p]] <- eval(parse(text=paste(panel.types[p],'_build(plots[[p]], p, DF, a)',sep='')))


  } else {

	stop(paste("unknown panel type -- '",panel.types[p],"'", sep=''))	# if no such function exists, lmPlot2 errors out 

  }

    # reset DF table (delete added tmp.data columns and what not
  DF <- DF.hold	

 }


##############################
##### construct the plot #####
##############################

lwidth <- pps
for(w in 1:length(all_atts(a,'panel.width'))) lwidth <- c(lwidth, all_atts(a,'panel.width')[w], pps)

  # sets layout according to specified widths in function arguments
lmLayout <- grid.layout(nrow = 1, ncol = length(lwidth), 
			widths = unit(lwidth, rep("null", length(lwidth))), 
			heights = unit(rep(1, length(lwidth)), rep("null", length(lwidth))))
plots$layout <- lmLayout
plots$plot.width <- plot.width
plots$plot.height <- plot.height


class(plots) <- "mm"

if(print.file=="no") print.file <- NULL
print(plots, name=print.file, res=print.res)


invisible(plots)


}   ###### END FUNCTION ######

