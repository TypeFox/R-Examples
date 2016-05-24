# routines for creating flow allocation matrices, classification rasters and hsu group
# summaries for the Dynamic TOPMODEL
require(raster)
#require(igraph)
#require(topmodel)
require(rgdal)

# return a matrix of indexes of<- lls immediately downslope of the given cell in the first column
# and correspoinding flow proportion allocated according to weighted midpoint slope in teh second
DownslopeFlowAllocations <- function(rast, cur,
                                     thresh=0)  # thresh is maximum slope to allow downslope flow
  # +ve value allows flow to go "uphill"
{
  # only consider raster cell that have a classification attached
  rast2 <- raster::setValues(rast, NA)
  rast2[cur]<-rast[cur]
  # ensure there are flow directons for all the elevation cells considered
  rast <- fill.sinks(rast, deg=0.1, silent=F)
  # adjacent cells
  # cur <- cur[which(!is.na(rast[]))]
  adj <- adjacent(rast,cur,directions=8, pairs=T)
  # select only directions with -ve (downslope) flow
  dz <- rast[adj[,2]] - rast[adj[,1]]
  good <- which(!is.na(dz) & dz <= thresh)
  adj <- cbind(adj[good, ], dz[good], NA)

  # divvie up flow direction by cell
  cells <- split(adj, as.factor(adj[,1]))
  ncells <- length(cells)
  #	pb <- txtProgressBar(max=length(ncells), title="flow allocations", style=3)
  #	on.exit(close(pb))
  adj <- lapply(cells,
                function(cell.dirs)
                {
                  # row index
                  #   setTxtProgressBar(pb, which(adj[,1]==names(cell.dirs))/ncells)
                  # rebuild the destination cell matrix
                  adj <- matrix(cell.dirs, ncol=4)
                  # sum
                  dz.cell <- adj[,3]  #adj[i.adj.cell,3]
                  if(all(dz.cell==0))
                  {
                    #browser()
                    # deal with flat areas by allocating equally in all directions
                    p <- 1/length(dz.cell)
                  }
                  else
                  {
                    # calculate weighted averages in each direction out of
                    p <- abs(dz.cell/sum(dz.cell))
                  }
                  adj[,4]<-p
                  return(adj)
                }
  )
  # co-erce back to a table and remove any nulls
  adj <- do.call(rbind, adj)
  # table: first col is source, second destination, third proportion of flow in that direction
  return(adj[,c(1:2,4)])
  # no downslope cells
  return(NULL)
}


# construct a nhru x nreach matrix, the (i,j)th entry gives the proportion of the baseflow out of
# the jth land HRU into the ith reach
build.chan.dist.matrix <- function(dem,
                                   drn,
                                   hru,       # raster of groupings
                                   chan.width=2)
{
    # build a raster identifying all reaches by ID
    reaches <- build.reach.raster(dem, drn, chan.width=chan.width)
    # use this determine how baseflow gets allocated between reaches
    w.chan <- build.flow.dist.matrix(dem, cm=hru, reaches=reaches)
    nreach <- nrow(drn)
    # just the land to channel transitions
    w.chan <- w.chan[-(1:nreach),-((nreach+1):ncol(w.chan))]
    return(normalise.rows(w.chan))
}

# create raster of channel locations
build.reach.raster <- function(dem,
         drn, nchan=nrow(drn),
         copy.to.mem=T,
         atb=NULL,
         atb.thresh=0.8,  # if drn not supplied then use this as threshold contributing area
         chan.width=5)
{
  # build the raster of river cells. Value gives the proportion of ech cell occupied by the channel
  if(is.null(drn) & !is.null(atb))
  {
    a <- upslope.area(dem, fill.sinks=T)
    message("Identifying channel(s) from TWI...")
    a.thresh <- max(atb[], na.rm=T)*atb.thresh
    # use the TWI to idenifty the channel as those areas exceeding the threshold
    reaches<- atb > a.thresh
    reaches[which(reaches[])] <- 1
    reaches[which(reaches[]==0)] <-NA
  }
  else if(!is.null(drn))
  {
    if(copy.to.mem)
    {
      dem <- rast.mem.copy(dem)  # to prevent disk thrashing
    }
    reaches <- dem
    reaches[] <- NA
    # one row per channel
    reach.cells <- lapply(extract(dem, drn, cellnumbers=T),
                          function(r)
                          {
                            r[,1]
                          })

    ilen <- sapply(reach.cells, length)
    # ensure each reach identified with at least cell, so apply in reverse order of length
    for(i in order(sapply(reach.cells, length), decreasing=T))
    {
      reaches[reach.cells[[i]]] <- i
    }
  }
  else
  {
    stop("Supply channels as shapefile or supply raster of TWI")
  }

  # calculate the cell proportions
  prop <- min(chan.width/xres(dem), 1)
  #	cellprops <- reaches*min(chan.width/xres(dem), 1)

  # add a layer with proportions of cell occuppied by channel. estimated by
  # proportion of cell size to channel width, probably close enough
  reaches <- addLayer(reaches, reaches*prop)

  names(reaches)=c("chan", "chanprops")

  return(reaches)
}

do.build.flow.matrix <- function(dem, hru, drn, chan.width=4, fact=2)
{
	message("aggregating....")
	dem <- aggregate(dem, fact)
	hru <- aggregate(hru, fact, fun=modal)
	message("getting reaches")
	reaches <- build_chans(dem, drn, chan.width=chan.width)
	message("building...")
	w <- build.flow.dist.matrix(dem, hru, reaches=reaches)
}

# Create a weighting matrix from a similar dem and raster of cell classifications
#                   raster built from spatial object representing channel. If supplied, it should be a spatialines
#                   data frame with a row for each  reach. In the weighting matrix HSU ids are
#                   created for each reach and inserted before the landscape
#                   HSUs. Landscape cells that contain part of the channel flow to that reach,whatever the local topography
build.flow.dist.matrix <- function(dem, cm,  # landscape units
                  drn=NULL,
                  reaches=NULL,
									ndp=3,
                  reach.outputs=T,   # if true row and cols added for transitions out of channel units.otherwise only tarnsitions into channel units from land considered
								  all=F)
{
  cm <- cm + 1000
#  dem <- rast.mem.copy(dem)
#  reaches <- rast.mem.copy(reaches)
#  cm <- rast.mem.copy(cm)
	# check	that rasters have same dims and resolution
#	compareRaster(dem, cm, res = TRUE)
	# note: these are the ids that indicate HSU classifications, which may contain
	# information on how the hsu was discretised. The transition matrix
	# assumes that they are in sequential order, which can be obtaied from the sort order
	# as teh discretisation retains the structure of the groupings in the order
	# save original ids
	if(all)
	{
		# include every cell, including those outside catchment
		hsuids <- 1:length(cm)
		cellnos <- hsuids
	}
	else
	{
		# just non-NA
		hsuids <- unique(cm[[1]], na.rm=T)  # 1:max(cm[], na.rm=T)}
		cellnos <- which(!is.na(cm[]))
	}
  #back up unit names excluding the channels
  land.ids <- hsuids-1000
	#reaches <- build.chans(dem, drn, reaches, chan.width=4)
	reachids <- setValues(cm, 0)
  cellprops <- setValues(cm, NA)

	# add in channel cells if
  if(!is.null(reaches)) # & !is.null(drn))
  {
    reachids <- reaches[[1]]
    # proportion of cells occupied by reach
    cellprops <- reaches[[2]]
    # reorder so that reach ids are sequential
    rids <- unique(reaches[[1]])

 #   reachids <- subs(reaches[[1]],
#                     data.frame(rids, order(rids)))
    # insert a hsu for every reach without duplicates. add prefix to prevent duplication
    # with land units

    hsuids <-  c(1:max(rids), hsuids)
  }

	# reclass raster so sequential ids are used
	cm <- raster::subs(cm, data.frame(hsuids, order(hsuids)))

	cat("Getting downslope flow weights...\n")
  start.tm <- Sys.time()
  down.all <- DownslopeFlowAllocations(dem, cellnos)
	cat(nrow(down.all), "directions processed in ", format(difftime(Sys.time(), start.tm), digits=2), "\n")
 # down.all <- down.all[which(!is.na(down.all[,2])),]

  if(length(down.all)==0)
  {
    warning("no flow paths identified")
    return(matrix(0, nrow=length(hsuids), ncol=length(hsuids)))
  }
  from <- down.all[,1]
  to <- down.all[,2]
	down.all <- cbind(down.all, cm[from], cm[to])

  # add in river transitions
	reach.cells <- which(reachids[]>0)
  ireach <- which(down.all[,2] %in% reach.cells)
  rtrans <- down.all[ireach, ]
  # proportion of flow from these cells into channel (rather than identified HSU transition)
  props <- cellprops[rtrans[,2]]
  # build extra rows two for each cell- river cell transition
  # split flow according to proportion of cell occupied by channel
	riv.props <- props*rtrans[,3]   # third col is proportion
  land.props <- (1-props)*rtrans[,3]
  # replace with updated land transitions and append river transitions
  down.all[ireach,3] <- land.props
  # create new rows. 1 is from cell, 2 dest cell. Replace col 3, prop, with calculated value
  # from HSU teh same, last col is destination HSU id - river ids come before HSU ids
  new.rows <- cbind(down.all[ireach, 1:2], riv.props, down.all[ireach,4], reachids[rtrans[,2]])
  down.all <- rbind(down.all, new.rows)
  # HSU and channel transistion table: from hsu, to hsu (or channel), flow proportion
#  trans <- data.frame(as.factor(down.all[,4]), as.factor(down.all[,5]), down.all[,3])

  trans <- data.frame(down.all[,4], down.all[,5], down.all[,3])
 trans[,1]  <- as.factor(hsuids[trans[,1]])
 trans[,2]  <- as.factor(hsuids[trans[,2]])
  names(trans)<-c("from", "to", "prop")
  trans <- trans[which(!is.na(trans[,2])),]

  # missing destinations ie units that are not linked to any others
  # cross tabulate i.e collate tranistion between pairs of elemenets into a table nhruxnhru in size
  cat("Cross tabulating inter-group transitions...\n")
  w <- xtabs(prop~from+to, data=trans)

  # add in any missing transitions
  while(length(which(!hsuids %in% colnames(w)))>0)
  {
    imissing <- which(!hsuids %in% colnames(w))[1]
    nm <- hsuids[imissing]
    #needs looking at
    w <- insert.col(w, imissing, 0, nm)

  }

  if(!is.null(reaches))
  {
    colnames(w) <- c(paste0("R", 1:max(rids)), land.ids)
  }
  else
  {
    colnames(w) <- land.ids
  }

  zero.rows <- which(rowSums(w)==0)
  if(length(zero.rows)>0)
  {
    # transitions into nowhere!

    warning(paste0("Nil row sum in flow matrix", paste0(names(w)[zero.rows], collapse=", ")))
   # w <- w[-zero.rows,]

  }
# w <- signif(w, ndp)
	# round to sensible no. dp and renormalise to rows add to 1

#

 	if(!is.null(reaches))
 	{
 	  nreach <- max(rids)  #might actually be fewer
     if(reach.outputs)
     {
      # insert an identity matrix for channel transistions
      # and route all channel flow via a time delay procedure and /or construct
      # inter channel flow transition matrix to route flows down channel
       rw <- matrix(0, ncol=ncol(w), nrow=nreach)
       rw[1:nreach, 1:nreach] <- identity.matrix(nreach)
       w<-rbind(rw, w)

       rownames(w)<- colnames(w)
     }
     else
     {
       # only the land to river transitions
       w <- w[,1:nreach]
     }
	}
  else
  {
    rownames(w) <- land.ids
  }


  w <- normalise.rows(w)

	return(w)
}


insert.col <- function(mat, i, val=NA, nm=NULL)
{
  if(length(val)==1){val<- rep(val, nrow(mat))}
  res <- cbind(mat[,1:min(i, ncol(mat))], val)
  colnames(res)[ncol(res)]<- nm
  if(i<ncol(mat))
  {
    res <- cbind(res, mat[,(i+1):ncol(mat)])
  }

  return(res)
}

normalise <- function(x)
{
	if(is.vector(x))
	{
		tot <- sum(abs(x), na.rm=T)

		return(ifelse(tot==0,x,x/tot))
	}
	return(x)
}




# follow the path downslope until reaching
get.paths.to.outlet <- function(pths, ipths=NULL, outlet)
{
  if(length(ipths)==0)
  {
    # cat("finished")
    return(NULL)
  }
  res <- ipths
  for(ipth in ipths)
  {
    inxt <- which(pths$NODE_B==pths[ipth,]$NODE_A)
    res <- c(res, get.upslope.paths(pths, inxt))
  }
  return(unique(res))
}

get.upslope.paths <- function(pths, ipths=NULL)
{
  if(length(ipths)==0)
  {
   # cat("finished")
    return(NULL)
  }
  res <- ipths
  for(ipth in ipths)
  {
    inxt <- which(pths$NODE_B==pths[ipth,]$NODE_A)
    res <- c(res, get.upslope.paths(pths, inxt))
  }
  return(unique(res))
}


# straight line distance. could scale up to get a nearer estimate
simple.dist.to.outlet <- function(dem, outlet, scale.fact=1.5)
{
	if(!is.null(outlet))
	{
			dists <- distanceFromPoints(dem, xyFromCell(dem, outlet))+ dem-dem
			return(dists*scale.fact)
	}
}

# determine cell in region of cell, if supplied, or lowest cell in DEM if not, that has
# the most upslope connectivity
locate.outlet <- function(dem, cell=NULL)
{
	bound <- raster::boundaries(dem)
	bound[bound==0] <- NA
	# only cells on edge
  outlet <- which.min((dem*bound)[])

 	return(outlet)

  adj <- adjacent(dem,cell,directions=8,pairs=F, include=T)

  ndest <- sapply(adj,
        function(x)
        {
          adj <- adjacent(dem,x,directions=8,pairs=F)
          upslope <- adj[which(dem[adj]>dem[x])]
          ndest <- length(upslope)
        }
  )

  outlet <- adj[which.max(ndest)]
  return(outlet)

}

# flow.lens.2 <- function(dem, cells=NULL, dest=NULL, samp=30)
# {
#   dem <- fill.sinks(dem, deg=0.5, silent=F)
#   dir <- terrain(dem, "flowdir")
#   if(is.null(cells))
#   {
#     cells <- which(!is.na(dem[]))
#   }
#   cells <- sample(cells, size=30)
#   pb <- txtProgressBar(max=length(cells), style=4)
#   i <<- 1
#   pths <- sapply(cells,
#          function(cell)
#         {
#           pth <- flowPath(dir, cell)
#           i <<- i + 1
#           setTxtProgressBar(pb, i)
#           return(pth)
#
#         })
#   if(!is.null(dest))
#   {
#     dests <- lapply(pths, function(pth) { pth[length(pth)]==dest})
#   }
#
#
# }

# MatchClass <- function(x, classes)
# {
# 	return(class(x)[1] %in% classes)
# }

# use the lengths and dem to obatin a 3-band grouping foet he catchment
# then apply the no. custs to teh corresponding layer. Combine to
# produce a final HSU classification
# if channels are supplied then remove the river cells from the calculation
combine.groupings <- function(dem,
                             	layers=list(),
														 	catch=NULL,
														 	chans=NULL,
							 								cuts=c(a=5),
              								thresh=2/100,   # threshold hsu area contrib
														 	equal.areas =F ) # if T then groups are strictly equal in plan area
{
	if(is.null(catch))
	{ # has to be at least one layer supplied
		layers <- delete.NULLs(layers)
  	catch<- raster::stack(layers)
	}

  if(!is.null(chans))
  {
    # remove river cells from calcs
    ichan <- which(chans[[1]][]>0)
    catch[ichan] <- NA
  }

  # select names in stack matching discretisation
  nms <- intersect(names(catch), names(cuts))

  if(length(nms)==0){stop("No layers found in catchment raster stack matching names of cuts to be applied")}

  # default is just one HSU representin entire catchment
	cm <- dem-dem+101

	ints <- list()
	#if(!all(names(cuts) %in% names(catch)))
	# apply hsu id for point by combining the cut value in each layer
	for(nm in nms)
	{
		for(id in raster::unique(cm))
		{
			idcells <- which(cm[]==id)
			n.cut <- cuts[[nm]]
			# cut cells already in this grouping according to the sublevel
			if(is.null(n.cut))
			{
				stop(paste(nm, " specified as cut variable but no corresponding layer found"))
			}
			if(n.cut>1)
			{
				laycutvals <- cut(catch[[nm]][idcells], n.cut, labels=F)


			}
			else{
        # one cuts only so return just when a nin-NA values existin
				laycutvals <- as.numeric(catch[[nm]][idcells]>=0)
			}
			# build new id from top level category plus sub level
			cm[idcells]<-10*cm[idcells] + laycutvals
		}
	}
  # need this to distinguish channel and land HSUs
  if(max(cm[],na.rm=T)<100)
  {cm<-cm+100}

  # merge smaller groups
  # source list
  cm <- merge.groups(cm, thresh)
  cm <- merge.groups(cm, ids=rev(unique(cm)), thresh)

  subs2 <- data.frame(cbind(unique(cm, na.rm=T), 100+order(unique(cm, na.rm=T))))
  cm <- subs(cm, subs2)


	names(cm)<-"HRU"
	# add in the dem and each discrteised layer
	cm <- addLayer(cm, catch)

	return(cm)
}

get_group_bounds <- function(cm)
{
    nms <- setdiff(names(cm), "HRU")
    hru <- cm[["HRU"]]
    # determine set bounds
    bnds <- NULL
    for(nm in nms)
    {
        cats <- split(cm[[nm]][], hru[])
        cats <- lapply(cats, range, na.rm=T)
        cats <- do.call(rbind, cats)
        colnames(cats)<-paste(nm, c("min", "max"), sep=".")
        bnds <- cbind(bnds, cats)
    }

    return(bnds)
}

merge.groups <- function(cm, thresh, ids=unique(cm))
{
 # ids <- unique(cm)
  # id mapping
  map <- data.frame()
  i <- 1
  ntot <- length(which(!is.na(cm[])))
  while(i <= length(ids))
  {
    id <- ids[i]
    tot <- length(which(cm[]==id))
    j<-1
    map <- rbind(map, c(id, id))
    while((i+j)<=length(ids) & (tot/ntot)<thresh)
    {
      idj <- ids[i+j]
      tot <- tot + length(which(cm[]==idj))
      map <- rbind(map, c(idj, id))
      j<-j+1
    }
    # keep collecting until reaching threshold
    i<- i+j
  }
  ids <- rev(ids)
  return(subs(cm, map))
}

# PlotGroupings <- function(cm, dem, drn, sel=extent(cm),legend=T, gridcells=F, ...)
# {
#
# 	groups <- build.hru.table(cm)
# 	ngroups <- nrow(groups)
# 	nchan <- length(which(groups[,1]<100))
# 	cm <- cm[[1]]
# 	# substitute values so group numbers are sequenetial
# 	cm <- subs(cm, groups)
# 	cm<-cm-nchan
# 	cm[which(cm[]<0)]<-0
# 	cols <-  terrain.colors(ngroups-nchan)
# 	#(rep("blue", nchan),
# 	PlotDemADrn(dem, drn, a=cm, sel=sel, col=c("lightblue", cols), legend=F, ...)
# 	if(legend)
# 	{
# 	legend(x="bottomleft", legend=groups[(nchan+1):nrow(groups),3],
# 		   ncol=2,bg="white",fill=cols, title="Eff distance to channel (m) / log(a) / slope angle ?")
# 	}
#
# 	if(gridcells)
# 	{
# 		HighlightCells(cm,cellsFromExtent(cm,sel))
# 	}
#
# }
#
#
# # draw arrows from a cell to downslope location with thicknesses proportional
# #to teh proportion of flow in each directiom
# ShowAllocation <- function(dem, cellno, maxwidth=5, cex.text=1, label=T)
# {
# 	down <- DownslopeFlowAllocations(dem, cellno)
#
# 	fromxy <- xyFromCell(dem, cellno)
#
# 	apply(down, MARGIN=1,
# 		  FUN=function(tocell)
# 		  {
# 		  	width <- max(1,tocell[2]*maxwidth) # < maxwidth
# 		  	toxy <- xyFromCell(dem, tocell[1])
# 		  	arrows(fromxy[1],fromxy[2],toxy[1],toxy[2],lwd=width,length=0.1)
# 		  	if(label)
# 		  	{
# 		  	mid <- (fromxy+toxy)/2
# 		  	text(mid[1], mid[2], round(tocell[2],1), cex=cex.text)}
# 		  }
# 	)
# }

build.hru.spatial <- function(disc, drn=disc$drn)
{
  if(!is.null(drn))
  {
    # convert drn to polygons
    drn <- rgeos::gBuffer(drn, width=max(round(disc$chan.width/2), 1), byid=T)
  }
#  drn <- SpatialPolygonsDataFrame(drn, data=data.frame())
  # create a spatial polygons dataframe object from the given
  # for the given discretisation
  ids <-  disc$groups[,"id"]
  ngroups <- length(ids)  # includes river

  hru.rast <- disc$hru.rast

  cat("Building HRU polygons...")
  hru.sp <- rasterToPolygons(disc, dissolve=T)
  polys <- hru.sp@polygons
  poly.ids <- unique(hru.rast)
    #     if(!is.null(drn))
    #     {
    #       drnbuff <- gBuffer(drn, width=max(round(chan.width/2+0.5), 1), byid=T)
    #       polys <- c(drnbuff@polygons, polys)
    #       poly.ids <- c(1:length(drnbuff@polygons), poly.ids)
    #     }
    # reconstruct SpatialPolygons object
  hru <- SpatialPolygons(polys, proj4string=hru.rast@crs)
  # check that group and polygon ids are equal
  try(if(!all(ids==poly.ids)){warning("error in groups and hru class ids")})
  hru <- maptools::unionSpatialPolygons(hru, IDs=poly.ids)

  dat <- data.frame(disc$groups[disc$groups[,"id"] %in% poly.ids,])
  row.names(dat) <- dat$id
  hru <- SpatialPolygonsDataFrame(hru, dat)

}

mean.na <- function(x){return(mean(x, na.rm=T))}

# shif the supplied raster so its BL corner is coincident woth the given origin
reset.origin <- function(rast, origin=c(x=0,y=0))
{
  shift(rast, origin[1]-extent(rast)@xmin, origin[2]-extent(rast)@ymin)

}

# group summary table with optional raster of cell proportions available for land
# optionally supply raster with proportions of cells occupied by land: defaults to 1
# if no channel

build.hru.table<- function(cm,
                       dem=NULL,
                       reaches=NULL,
                       cellareas=cm-cm+1,  # props occupied by land 0-1
                       catch=NULL)
{
	if(is.null(catch))
	{
		catch=cm[[2:nlayers(cm)]]
		# assume that catchment discretistaion info is passed via the other layer of
		# the multi-band rasteres
		cm <- cm[[1]]
	}
    a.atb <- NULL
    if(all(c("a", "atb") %in% names(catch)))
    {
	    a.atb <- catch[[c("a", "atb")]]
    }
	ids <- unique(cm[[1]])
 	# handle outside the class matrix
    cellareas[which(is.na(cellareas[])&!is.na(cm[]))]<-1
	# maximum plan area of land within each cell
  maxCellArea <- xres(cm)*yres(cm)

	cat("Building areas...")
	areas<- sapply(ids,
                 function(id)
                   {
                    cm.props <- cellareas[]*(cm==id)[]
             #       browser()
                    sum(cm.props[], na.rm=T)* maxCellArea
                 }
	)

	# add in river reaches, removing the area covered by channel from the
	# corresponding groups
	if(!is.null(reaches))
	{
		# determine land areas occupied by each channel and insert ids and areas
		# at head of group table
		chans <- zonal((1-cellareas)*maxCellArea, reaches, fun=sum) # sum the areas of cells classed by channel id multiplied by proportion ocuupied for edach
		# two colums: id, area (nas removed)
		areas <- c(chans[,2], areas)
		ids <- c(chans[,1], ids)
# 		ids <- c(rids, ids)
	}

  # total catchment area (includes reaches)
  totArea <- sum(areas)

	# reordering so that river reaches appear first. add 100 to distinguish reaches from land hsus
#	orders <- c(100+((nchan+1):ngroups), 1:nchan)

	# build table, first two columns will be used to renumber the
	groups <- data.frame("id"=ids,
	                     "tag"=ids,
                       "chan.no"=NA,
						 "order"=1:length(ids),
		#				 "breaks"=ints,
						 "area_pc"=round(100*(areas/totArea),2),
						 "area"=round(areas))

#	groups <- cbind(groups, "atb.bar"=0)
  groups <- add.upslope.areas(groups, dem, cm, a.atb=a.atb,
  													area_pcs=cellareas)
	groups$atb.bar <- round(groups$atb.bar, 2)
	groups <- cbind(groups, "gauge.id"=1)
  groups <- cbind(groups, "catch.id"=1)
    add.par <- def.hsu.par()
    nms <- setdiff(names(add.par), names(groups))
    add.par <- add.par[nms]
  pars <- data.frame(matrix(rep(add.par, nrow(groups)), byrow=T, nrow=nrow(groups)))
  colnames(pars) <- nms

  groups<- cbind(groups, pars)
  groups <- apply(groups, MARGIN=2, FUN=function(x){unlist(x)})

    #row.names(groups) <- groups[,1]
	return(data.frame(groups, 	row.names=groups[,1]))
}

# locate nearest aws to each group and return index
add.nearest.gauge <- function(groups, hru.sp, drn, gauges)
{
  ids <- groups$id
#  groups <- sort(groups$hru.sp)

  for(i.group in 1:nrow(groups))
  {
    id <- groups[i.group,]$id
    hru.geom
    # locate the geometry associated with the group
    hru.geom <- hru.sp[which(hru.sp$HRU==id),]
    cent <- rgeos::gCentroid(hru.geom)
    # locate current max
    dist <- rgeos::gDistance(cent, gauges[i.gauge,])
    for(i.gauge in 1:nrow(gauges))
    {
      if(rgeos::gDistance(cent, gauges[i.gauge,])<dist)
      {
        groups[i.group,]$gauge.id  <- i.gauge
      }
    }
  }
  hru.sp@data <- groups
  return(hru.sp)
}

# add total upslope area and mean topographic index to groups info
add.upslope.areas <- function(groups, dem, class.m, a.atb=NULL,
														area_pcs=round(dem/dem))   # optional ratser of cell areas occupied by land
{
  if(!(is.null(dem) | is.null(class.m)))
  {
  	if(is.null(a.atb))
  	{
    	# mean ln(a/tan(b))
    	a.atb <- upslope.area(dem, atb=T)
  	}
# uncomment to remove cells containg channel from analysis
      area_pcs[area_pcs<1]<- NA
    # deal with cells partly ocuppied by river channel
    atb.adj <- a.atb[["atb"]]+log(area_pcs)
    # zonal statistics - mean is default. Values adjusted for any cells containing
    # channel
    atb <- zonal(atb.adj, class.m, "mean")  # deals with nas

    # specific discharge per unit recharge [-] assuming steady state.
    # a measure of the relative yield of the area?
  #  sigma.a <- zonal(raster::setValues(dem,a.atb$area), cm, "mean")  # this must have the same first row
    for(row in 1:nrow(atb))
    {
      id <- atb[row,1]
      indx <- which(groups$"id"==id)
      if(length(indx)==0)
      {
        warning(paste("index ", id, "from class matrix not found in groups table"))
      }
      else
      {
      #  cat(id, "\t", indx, "\n")
      # this assummes that every id in the raster has a corresponding id in the table....
      groups[indx,"atb.bar"]<- atb[row,2]
      # q over r
     # groups[indx,"sigma.a"]<- sigma.a[row,2]/groups[indx,]$area
      }
    }
  }
  return(groups)
}

mem.copy <- function(rast)
{
  if(raster::inMemory(rast)){return(rast)}
  #  rast.copy <- stack(rast)
  layers <-
    lapply(names(rast),
           function(ln)
           {
             layer <- rast[[ln]]
             raster::setValues(layer, getValues(layer))
           }
    )
  if(length(layers)>1)
  {
    return(stack(layers))
  }
  else
  {
    return(layers[[1]])
  }
}
#############################################################
# raster of reach locations and cell proportion occupied
#############################################################
build.reach.raster <- function(dem, drn, nchan=nrow(drn),
															 copy.to.mem=T,
															 atb=NULL,
															 atb.thresh=0.8,  # if drn not supplied then use this as threshold contributing area
                               chan.width=1)
{
	# build the raster of river cells. Value gives the proportion of ech cell occupied by the channel
	if(is.null(drn) & !is.null(atb))
	{
		a <- upslope.area(dem, fill.sinks=T)
		message("Identifying channel(s) from TWI...")
		a.thresh <- max(atb[], na.rm=T)*atb.thresh
		# use the TWI to idenifty the channel as those areas exceeding the threshold
		reaches<- atb > a.thresh
		reaches[which(reaches[])] <- 1
		reaches[which(reaches[]==0)] <-NA
	}
	else if(!is.null(drn))
	{
		if(copy.to.mem)
		{
			dem <- mem.copy(dem)  # to prevent disk thrashing
		}
		reaches <- dem
		reaches[] <- NA
		# one row per channel
		reach.cells <- lapply(extract(dem, drn, cellnumbers=T),
													function(r)
													{
														r[,1]
													})

		ilen <- sapply(reach.cells, length)
		# ensure each reach identified with at least cell, so apply in reverse order of length
		for(i in order(sapply(reach.cells, length), decreasing=T))
		{
			reaches[reach.cells[[i]]] <- i
		}
	}
	else
	{
		stop("Supply channels as shapefile or supply raster of TWI")
	}

	# calculate the cell proportions
	prop <- min(chan.width/xres(dem), 1)
#	cellprops <- reaches*min(chan.width/xres(dem), 1)

	# add a layer with proportions of cell occuppied by channel. estimated by
	# proportion of cell size to channel width, probably close enough
	reaches <- addLayer(reaches, (reaches>0)*prop)

	names(reaches)=c("chan", "chanprop")

	return(reaches)
}


get.def.thresh <- function(area.thresh=1)
{
  return(area.thresh)
}

get.defs <- function()
{
  return(list(chan.width=2, area.thresh=1))
}


# create a new discretisation and add to project's list. if it exists already then just add
# to list of discretisations, unless rebuild is true in which case rebuid the components
# add.disc <- function(proj, cuts=NULL, i.disc=0, rebuild=F, chan.width=2,
# 										 ...)
# {
# 	# proj <- merge.lists(proj, list(...))
#
# 	# build names from extra parameters
# 	cuts <- merge.lists(cuts, list(...))
#
# 	# name wil be infered from cuts
# 	dn <- disc.dir.name(cuts, dn=proj$disc.dir)
# 	catch <- NULL
# 	try(catch <- build.proj.catch(proj))
# 	# rebuild <- file.exists(dn)
# 	if(rebuild){message(paste("Rebuilding files in ", dn))}
# 	disc <- NULL
# 	disc <-disc.from.dir(dem=proj$dem,
# 													 drn=proj$drn,
# 													 reaches=proj$reaches,
# 													 routing=proj$routing,
# 													 catch=catch,
# 													 dn=dn,
# 													 cuts=cuts,
# 													 rebuild=rebuild,
# 													 chan.width=chan.width,
# 													 agg=agg,
# 												#	 area.thresh=proj$area.thresh,
# 												...)
#
# 	if(!is.null(disc))
# 	{
# 		if(is.null(i.disc) | i.disc <=0 | i.disc > length(proj$disc))
# 		{
# 			# return new discretisation
# 			return(disc)
# 		}
# 		# same cuts as any other?
# 		proj$disc[[i.disc]]<- disc
# 		return(proj$disc[[i.disc]])
# 	}
# 	else
# 	{
# 		warning("Adding discretisation ", paste(cuts, collapse=","), "failed")
# 	}
# }




#*******************************************************************************
# TOPMODEL helper routines
#*******************************************************************************
# discretisation of catchment according to topographic index
disc.topidx <- function(dem, riv, nbreaks)
{
  atb<-upslope.area(dem, atb=T)[["atb"]]
  if(!is.null(riv))
  {
    atb[which(riv[]>0)]<-NA  # ignore river cells
  }

  atb.hist <- hist(atb, breaks=nbreaks, freq=F)
  tot <- sum(atb.hist$counts)
  # second col is for proportion occupied - final entry is zero for count of cells > final break
  topidx.frame <- cbind(atb=atb.hist$breaks, prop=c(atb.hist$counts/tot, 0))
  return(topidx.frame)
}


# build an approximate routing table from the dem
#"Flow is routed through a delay function which represents the time spent in the
# channel system. The parameter delay is used for this. Delay is a matrix with 2
# columns. The first column gives the cumulative relative area. The second column
# gives the average distance towards the outlet (m)."
get.routing.table <- function(dem, nreach=5)
{
  dem <- raster::setValues(dem, topmodel::sinkfill(as.matrix(dem), degree=0.1, res=10))

  fl <- raster::setValues(dem, topmodel::flowlength(as.matrix(dem)))*xres(dem)

  routing.hist <- hist(fl, breaks=nreach)
  n <- sum(routing.hist$counts)
  # build the table@ first column is disatnce to outlet, second the cumlative area
  routing <- cbind(prop=cumsum(routing.hist$counts)/n, d=routing.hist$mids)
  return(routing)
}
#*******************************************************************************

# construct a name for a directory in whci to put files for catchment discretisation
disc.dir.name <- function(cuts, dn="", area.thresh=0)
{
  # convert to %age
#  if(area.thresh<1){area.thresh <- area.thresh*100}
  return(file.path(dn, paste0(names(cuts),"=", cuts, collapse=",")))
}


is.stack <- function(obj)
{
  if(is.null(obj)){return(FALSE)}
  if(inherits(obj, "RasterStack") | inherits(obj, "RasterLayer") | is(obj, "RasterBrick")){return(TRUE)}
  return(FALSE)
}

# wrapper to create all input for Dynamic TOPMODEL run given DEM and a set of
# breaks on eff dist to channel, uplsloploe area and slope
# also builds information for "normal" TOPMODEL run

aggregate.null <- function(rast, fact, fun=mean)
{
  fact <- round(fact)
  if(is.stack(rast) & fact > 1)
  {
    rast <- aggregate(rast, fact, fun=fun)
  }
  return(rast)
}

get.flow.distribution.matrix <- function(dem, cm, reaches, sf=3)
{
	agg <- max(1,round((ncell(dem)/1e6)))
	if(agg>1){
		message(paste0("Aggregating dem by factor of ", agg))
	}
	dem <- aggregate.null(dem, agg)
	cm <- aggregate.null(cm, agg, fun=modal)
	reaches <- aggregate.null(reaches, agg)

	message("Creating flow transistion matrix....")
  # flow transistion matrix including transfers to river
	w <- build.flow.dist.matrix(cm=cm, dem=dem, reaches=reaches)
	# always require some kind of channel identification raster even if constructed from DEM
	message("Creating channel routing matrix....")
	# construct flow connectivity graph using the reach identification raster previously loaded
	# or built, using reach ids as the "classification" and flow weightings from dem.
	adj <- build.flow.dist.matrix(cm=reaches[[1]], dem=dem, reaches=NULL)

	# clear the weights for all river hsus, then add the adjacency matrix
	# overwrite channel entries of weighting matrix with river connections
	w<- as.matrix(w)
  ichan <- 1:nrow(adj)
	w[ichan,]<-0

	w[ichan, ichan] <-  as.matrix(adj)
	w <- signif(w, sf)
  row.names(w)[ichan]<- paste0("R", row.names(adj))
	ids <- row.names(w)

	rownames(w)<- ids
	colnames(w)<- ids
	return(w)
}

#
# create.reach.info <-  function(dem, w, reaches)
# {
#   # destinations after 100 cell-cell iterations
#   w100<-matrix.power(w, 100)
#
#   # consider just land - river transitions
#   w100<- as.matrix(w100[(nchan+1):nrow(w),1:nchan])
#
#   # normalise to a probability distribution - asuume all flow enters a reach
#   props <- round(NormaliseVector(colSums(w100)), 3)
#
#   # flow lengths to outlet. doesn't work very well without first filling sinks
#   flowlens <- raster::setValues(dem, flowlength(as.matrix(dem))*xres(dem))
#
# 	# in-channel distances determined by zonal mean of flow lens for dem cells containing the channel(s)
# 	chanlens <-zonal(flowlens, reaches[[1]])
#
#   # construct reach info table. Compatible with TOPMODEL if props are cumulative
#   reach.info <- data.frame("dist"=round(chanlens[,2]), "prop"=props)
#   return(reach.info)
# }
#
#




# not needed: Laplacian filter for raster focal function
#fact <-1/sqrt(2)
#downslope.contour.filter <- matrix(c(fact, 1, fact, 1,1, 1,fact,1,fact), nrow=3, byrow=T)

FindNext <- function(x)
{
 # browser()
  cur <- x[5]
  if(is.na(cur))
  {
    return(NA)
  }
  # weight straight directions more heavily
  dz <- x[5]-x
  return(which.max(dz))
}


