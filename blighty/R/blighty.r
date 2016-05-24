# the blighty package bydata David Lucy - plots the coastline and some 
# internal features of the British Isles
# now does Ireland in V3
# Copyleft - David Lucy January 2006

blighty <- function(
		    place="set.British.Isles",
		    set="TRUE", 		# standard set or vector of objects
		    grid=FALSE,
		    xlimits, 
		    ylimits,
		    xpadding=0,
		    ypadding=0,
		    parcol=par("bg"),		# primary area colour
		    parbor=par("fg"),		# primary area border colour - "transparent" for invisible
		    parwdh=1,			# primary area border width
		    sarcol=par("bg"),		# secondary area colour
		    sarbor=par("fg"),		# secondary area border colour - "transparent" for invisible
		    sarwdh=1,			# secondary area border width
		    parang=NULL,		# angle of lines for primary areas
		    parden=NULL,		# lines per inch for primary areas
		    sarang=NULL,		# angle of lines for secondary areas
		    sarden=NULL,		# lines per inch for secondary areas
		    tlncol=par("fg"),		# colour of lines for non-area objects
		    tlnwdh=1,			# width of lines for non-area objects
		    grdcol=par("fg"),		# colour of grid lines
		    grdwdh=1)			# width of grid lines

{

# means the set specified is one of the standard sets of objects
# thus the objectnames can be read straight in from the set file
# and we don't have to test the objects to see whether they
# actually exist or not
	if(set == TRUE)
		{
		cat("\nPlotting ", place, " be patient ...\n")

		# get the vector of objectnames from the set file
		data(list = place)
		objectnames <- as.vector(get(place)[,1])
		non.existant.objects <- NULL

		# assign the number of objects in the set
		noobjects <- length(objectnames)

		# read in the data from the object files
		data(list = objectnames)

		cat("Data loaded ...\n")
		}


# means the set specified is not one of the standard sets of objects
# we now need to see whether the objects exist or not as place will be a
# vector of objects rather than a string representing a set of vectors of 
# objects
# we assume here that the object being refered to by place is an array of strings
# with no col.name
	if(set == FALSE)
		{
		cat("\nPlotting ", place, " be patient ...\n")
		objectnames <- NULL
		non.existant.objects <- NULL

		# important to assign the place to an actual vector
		# of names
		tmp.objectnames <- get(place)


		# make sure that objects appear one only
		tmp.objectnames <- unique(as.character(tmp.objectnames))
		tmp.noobjects <- length(tmp.objectnames)


		# index through the place names and see whether
		# the objects exist for them in the data
		for(ctr in 1:tmp.noobjects)
			{
			does.object.exist <- suppressWarnings(exists(data(list = tmp.objectnames[ctr])))

			if(does.object.exist == TRUE)
				{objectnames[length(objectnames) + 1] <- tmp.objectnames[ctr]}

			if(does.object.exist == FALSE)
				{non.existant.objects[length(non.existant.objects) + 1] <- tmp.objectnames[ctr]}
			}

		noobjects <- length(objectnames)
		cat("Data loaded ...\n")
		}


# make sure that there is at least one valid object to plot - if not exit
if(length(objectnames) < 1){stop("No objects to plot - try another set or vector of placenames")}

#print(objectnames)
#stop("TEST STOP")

	# calculate x limits if not specified by the user
	# these are calculated from the map objects themselves
	if(missing(xlimits))
		{
		xlimits <- vector(mode="numeric", length=2)

			for(ctr in 1:noobjects)
				{
					if(ctr == 1)
						{
						len <- length(get(objectnames[1])$x)
						xlimits[1] <- min(get(objectnames[1])$x[2:len])
						xlimits[2] <- max(get(objectnames[1])$x[2:len])
						}

					if(ctr > 1)
						{
						len <- length(get(objectnames[ctr])$x)
						if(xlimits[1] > min(get(objectnames[ctr])$x[2:len])){xlimits[1] <- min(get(objectnames[ctr])$x[2:len])}
						if(xlimits[2] < max(get(objectnames[ctr])$x[2:len])){xlimits[2] <- max(get(objectnames[ctr])$x[2:len])}
						}
				}
		}


	# calculate x limits if not specified by the user
	# these are calculated from the map objects themselves
	if(missing(ylimits))
		{
		ylimits <- vector(mode="numeric", length=2)

			for(ctr in 1:noobjects)
				{
					if(ctr == 1)
						{
						len <- length(get(objectnames[1])$y)
						ylimits[1] <- min(get(objectnames[1])$y[2:len])
						ylimits[2] <- max(get(objectnames[1])$y[2:len])
						}

					if(ctr > 1)
						{
						len <- length(get(objectnames[ctr])$y)
						if(ylimits[1] > min(get(objectnames[ctr])$y[2:len])){ylimits[1] <- min(get(objectnames[ctr])$y[2:len])}
						if(ylimits[2] < max(get(objectnames[ctr])$y[2:len])){ylimits[2] <- max(get(objectnames[ctr])$y[2:len])}
						}
				}
		}


# test the plotting limits anything less than 200km in either direction
# leads to poor maps as the point resolution isn't that great
if(max(xlimits) - min(xlimits) < 200){cat("Less than 200km East-West - poor map resolution - proceeding\n")}
if(max(ylimits) - min(ylimits) < 200){cat("Less than 200km North-South - poor map resolution - proceeding\n")}

# add some padding around the figure if padding is specified 
# stops the figure from being hard up against the frame
xlimits[1] <- xlimits[1] - xpadding
xlimits[2] <- xlimits[2] + xpadding
ylimits[1] <- ylimits[1] - ypadding
ylimits[2] <- ylimits[2] + ypadding

# calculate suitable plot limits to keep the map with a more or less
# truly square grid
lims <- sqlimits(xlimits, ylimits)

cat("Correct aspect ratio calculated ...\n")

# setup the plot extremes
par(usr = c(lims$xlims[1], lims$xlims[2], lims$ylims[1], lims$ylims[2]))

# add to the plot all requested objects
#for(ctr in 1:noobjects){points(get(objectnames[ctr]), type="l")}
	for(ctr in 1:noobjects)
		{
		# 1 - area objects such as islands - landmasses etc
		if(get(objectnames[ctr])$x[1] == 1)
			{
			len <- length(get(objectnames[ctr])$x)
			polygon(get(objectnames[ctr])$x[2:len],
				get(objectnames[ctr])$y[2:len],
				angle=parang,
				density=parden,
				col=parcol,
				border=parbor,
				lwd=parwdh)
			}
		# 2 - area objects which sit on top of other area objects - things such as lakes
		if(get(objectnames[ctr])$x[1] == 2)
			{
			len <- length(get(objectnames[ctr])$x)
			polygon(get(objectnames[ctr])$x[2:len],
				get(objectnames[ctr])$y[2:len],
				angle=sarang,
				density=sarden,
				col=sarcol,
				border=sarbor,
				lwd=sarwdh)
			}
		# 3 - linear objects such as rivers - roads - railway lines
		if(get(objectnames[ctr])$x[1] == 3)
			{
			len <- length(get(objectnames[ctr])$x)
			points(get(objectnames[ctr])$x[2:len],
			       get(objectnames[ctr])$y[2:len],
			       type="l",
			       col=tlncol,
			       lwd=tlnwdh)
			}
		}

# if gridding has been requested put on the OS coordinates and gridlines and box
if(grid == "TRUE")
	{
	cat("Requested gridding applied ...\n")
	box(col=grdcol, lwd=grdwdh)
	xmarks <- pretty(lims$xlims)
	ymarks <- pretty(lims$ylims)

	axis(1, at=xmarks, labels=TRUE)
	abline(v=xmarks, col=grdcol, lwd=grdwdh)
	axis(2, at=ymarks, labels=TRUE)
	abline(h=ymarks, col=grdcol, lwd=grdwdh)
	}

# send the information to the global environment
blighty.mapinfo <- list(lims$xlims, lims$ylims, objectnames, non.existant.objects)
names(blighty.mapinfo) <- c("xlims", "ylims", "objects.used", "objects.not.found")
assign("blighty.mapinfo", blighty.mapinfo, envir=.GlobalEnv)

# do a bit of tidying up
rm(list = objectnames, envir = .GlobalEnv)

cat("Map complete ...\n\n")
}
