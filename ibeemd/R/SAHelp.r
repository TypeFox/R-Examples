# Created by Mao-Gui Hu (humg@lreis.ac.cn), 2013.10.
require(sp)
require(rgeos)
require(deldir)
require(spdep)

ClipPolygonsDf <- function(spInputPolys, spClipBnd)
{
	spDiff <- SpatialPolygonsDataFrame(
			gIntersection(spInputPolys[1,], spClipBnd), 
			spInputPolys[1,]@data, match.ID=FALSE)			
	for(i in 2:length(spInputPolys))
	{	
		gd <- gIntersection(spInputPolys[i,], spClipBnd)
		gd <- SpatialPolygonsDataFrame(gd, spInputPolys@data[i,], match.ID=FALSE)
		gd@polygons[[1]]@ID <- as.character(i)
		spDiff <- rbind(spDiff, gd)
	}
	return(spDiff)
}

# Generate voronoi SpatialPolygonsDataFrame from SpatialPointsDataFrame
CreateVoronoi <- function(spPointsDF, bDelaunay=FALSE, bClipVor=FALSE, spClipBnd=NULL)
{
	crds <- coordinates(spPointsDF)
	z <- deldir(crds[,1], crds[,2])
	if(!bDelaunay)
		w <- tile.list(z)
	else
		w <- triang.list(z)
	polys <- vector(mode='list', length=length(w))
	for (i in seq(along=polys)) 
	{
		pcrds <- cbind(w[[i]]$x, w[[i]]$y)
		pcrds <- rbind(pcrds, pcrds[1,])
		polys[[i]] <- Polygons(list(Polygon(pcrds)), ID=as.character(i))
	}
	SP <- SpatialPolygons(polys)
	if(!bDelaunay)
	{
		dd <- spPointsDF@data
		row.names(dd) <- sapply(slot(SP, 'polygons'), function(x) slot(x, 'ID'))
		voronoi <- SpatialPolygonsDataFrame(SP, data=dd)
		
		if(bClipVor)
		{
			if(is.null(spClipBnd))
			{
				#gu <- GetUnionSPolyDf(CreateVoronoi(spPointsDF, bDelaunay=TRUE))
				hpts <- chull(coordinates(spPointsDF))
				hpts <- c(hpts, hpts[1])
				gu <- gBuffer(SpatialPolygons(list(Polygons(list(Polygon(coordinates(spPointsDF)[hpts,])), ID="0"))),
						width=mean(sort(spDists(coordinates(spPointsDF)))[1:(length(spPointsDF)*4)]))
				spClipBnd <- SpatialPolygonsDataFrame(gu, data=data.frame(ID=0), match.ID=FALSE)
			}
			if(length(spClipBnd) > 1)
			{
				spClipBnd <- GetUnionSPolyDf(spClipBnd)
			}
			voronoi <- ClipPolygonsDf(voronoi, spClipBnd)
		}
	}
	else
	{
		voronoi <- SpatialPolygonsDataFrame(SP, data=data.frame(ID=1:length(SP)))
	}
	
	return(voronoi)
}

# Convert SpatialPolygonsDataFrame into SpatialLinesDataFrame
SPolygonsDfToSLinesDf <- function(spPolysDf)
{
	PolygonsToLines <- function(polys)
	{
		lstLine <- list()
		for(i in 1:length(polys@Polygons))
		{
			p <- polys@Polygons[[i]]
			curLine <- Line(p@coords)
			attr(curLine, "ringDir") <- p@ringDir
			lstLine[[i]] <- curLine
		}
		return(Lines(lstLine, polys@ID))
	}

	lstLines <- list()
	for(i in 1:length(spPolysDf))
	{
		p <- spPolysDf@polygons[[i]]
		lstLines[[i]] <- PolygonsToLines(p)
	}
	spLines <- SpatialLines(lstLines, spPolysDf@proj4string)
	spLinesDf <- SpatialLinesDataFrame(spLines, spPolysDf@data, match.ID=FALSE)
	return(spLinesDf)
}

# Convert SpatialLinesDataFrame into SpatialPolygonsDataFrame
SLinesDfToSPolygonsDf <- function(spLinesDf)
{
	LinesToPolygons <- function(lines0)
	{
		lstPoly <- list()
		for(i in 1:length(lines0@Lines))
		{
			l0 <- lines0@Lines[[i]]			
			curPoly <- Polygon(l0@coords)
			lstPoly[[i]] <- curPoly
		}
		return(Polygons(lstPoly, lines0@ID))
	}
	
	lstPolys <- list()
	for(i in 1:length(spLinesDf))
	{
		lines0 <- spLinesDf@lines[[i]]
		lstPolys[[i]] <- LinesToPolygons(lines0)
	}
	spPolys <- SpatialPolygons(lstPolys, proj4string=spLinesDf@proj4string)
	spPolysDf <- SpatialPolygonsDataFrame(spPolys, spLinesDf@data, match.ID=FALSE)
	return(spPolysDf)
}

# Convert SpatialPolygonsDataFrame into SpatialPointsDataFrame (centroid)
SPolygonsDfToSPointsDf <- function(spPolysDf)
{
	PolygonsToPoints <- function(polys)
	{
		ctr <- gCentroid(polys[1,])
		for(i in 2:length(polys))
		{
			ctr <- rbind(ctr, gCentroid(polys[i,]))
		}

		return(ctr)
	}
	
	spPoints <- PolygonsToPoints(spPolysDf)
	spPointsDF <- SpatialPointsDataFrame(spPoints, spPolysDf@data)
	return(spPointsDF)
}

# Union all polygons into a big area polygon
GetUnionSPolyDf <- function(spPolysDf, tolArea=NA, bLine=FALSE)
{
	MinArea <- function(spPolysDf)
	{
		minarea <- Inf
		for(i in 1:length(spPolysDf))
		{
			carea <- spPolysDf@polygons[[i]]@area
			if(minarea > carea) minarea <- carea
		}
		return(minarea)
	}
	
	if(is.na(tolArea)) tolArea <- MinArea(spPolysDf)/10
	
	# union all polygons into a big area polygon
	gu <- gUnaryUnion(spPolysDf)
	gu <- SpatialPolygonsDataFrame(gu, data.frame(id=1), F)
	
	# boundary of the polygon
	spLinesDf <- SPolygonsDfToSLinesDf(gu)
	
	# delete sliver polygons
	validIndx <- c()
	for(i in 1:length(spLinesDf@lines[[1]]@Lines))
	{
		crds <- coordinates(spLinesDf@lines[[1]]@Lines[[i]])
		p <- Polygons(list(Polygon(crds)), ID=i)
		if(p@area > tolArea) validIndx <- c(validIndx, i)
	}
	spLinesDf@lines[[1]]@Lines <- spLinesDf@lines[[1]]@Lines[validIndx]
	
	if(bLine)
		return(spLinesDf)
	else
		return(SLinesDfToSPolygonsDf(spLinesDf))
}

# Remove holes of a SpatialLinesDataFrame
RemoveHole <- function(spLinesDf)
{
	j <- 1
	lstLines <- list()
	for(i in 1:length(spLinesDf@lines[[1]]@Lines))
	{
		if(attr(spLinesDf@lines[[1]]@Lines[[i]], 'ringDir') == 1)
		{
			lstLines[[j]] <- spLinesDf@lines[[1]]@Lines[[i]]
			j <- j+1
		}
	}
	return(SpatialLinesDataFrame(
			SpatialLines(list(Lines(lstLines,ID=0))), data.frame(id=1), FALSE)
			)
}	

# Remove sliver polygons
RemoveSliver <- function(spPolysDf, tolArea=NA)
{
	MinArea <- function(spPolysDf)
	{
		minarea <- Inf
		for(i in 1:length(spPolysDf))
		{
			carea <- spPolysDf@polygons[[i]]@area
			if(minarea > carea) minarea <- carea
		}
		return(minarea)
	}
	
	Sliver <- function(polys)
	{
		delIndx <- c()
		for(i in 1:length(polys@Polygons))
		{
			p <- polys@Polygons[[i]]
			if(p@ringDir == 1 && p@area < tolArea)
				delIndx <- c(delIndx, i)
		}
		if(!is.null(delIndx))
		{
			polys@Polygons <- polys@Polygons[-delIndx]
			if(length(polys@Polygons) == 0) polys <- NULL			
			return(list(hasSliver=TRUE, polys=polys))
		}
		else
		{
			return(list(hasSliver=FALSE, polys=NULL))
		}
	}
	
	if(is.na(tolArea)) tolArea <- MinArea(spPolysDf)/10
	
	delIndx <- c()
	for(i in 1:length(spPolysDf))
	{
		p <- Sliver(spPolysDf@polygons[[i]])
		if(p$hasSliver)
		{
			if(is.null(p$polys))
				delIndx <- c(delIndx, i)
			else
				spPolysDf@polygons[[i]] <- p$polys
		}
	}
	if(!is.null(delIndx)) spPolysDf <- spPolysDf[-delIndx,]
	
	return(spPolysDf)
}

# Are polygons located at the border
LocatedAtBoundary <- function(spPolysDf, bRemoveHole=FALSE)
{
	cat("Boundary polygons ...\n")
	
	spLinesDf <- GetUnionSPolyDf(spPolysDf, bLine=TRUE)
	if(bRemoveHole) spLinesDf <- RemoveHole(spLinesDf)
	
	# whether each polygon intersect with the area boundary
	bBoundary <- rep(FALSE, length(spPolysDf))
	pb <- txtProgressBar(min=1, max=length(spPolysDf)-1, width=50, style=3)
	for(i in 1:length(spPolysDf))
	{
		setTxtProgressBar(pb, i)
		bBoundary[i] <- gIntersects(spPolysDf[i,], spLinesDf)
	}
	close(pb)
	
	return(bBoundary)
}

# First-order neighbors' index for SpatialPolygonsDataFrame
NeighborPolys <- function(spPolysDf)
{
	cat("Neighbor polygons ...\n")

	lstNeighbor <- poly2nb(spPolysDf, snap=1e-3)
	class(lstNeighbor) <- "list"
	attributes(lstNeighbor) <- NULL
	return(lstNeighbor)
}

# Find high order neighbours
FindNeighbours <- function(FirstNbList, nbOrder=2)
{
	if(nbOrder == 1) return(FirstNbList)
	
	highOrderNbList <- list()
	highOrderNbList[[1]] <- FirstNbList
	for(i in 2:nbOrder)
	{
		highOrderNbList[[i]] <- list()
		for(j in 1:length(FirstNbList))
		{
			if(length(highOrderNbList[[i-1]][[j]]) == 1 && 
				highOrderNbList[[i-1]][[j]] == 0)
				next
			
			nbs <- c()
			for(k in highOrderNbList[[i-1]][[j]])
			{
				nbIndx <- FirstNbList[[k]]
				# not occurred in lower neighbours
				for(l in 1:(i-1))
				{
					nbIndx <- setdiff(nbIndx, highOrderNbList[[l]][[j]])
				}
				nbs <- c(nbs, nbIndx)
			}
			highOrderNbList[[i]][[j]] <- sort(setdiff(unique(nbs),j))
		}
	}
	
	highOrderNbList$total <- list()
	for(i in 1:length(FirstNbList))
	{
		indx <- c()
		for(j in 1:nbOrder)
		{
			indx <- c(indx, highOrderNbList[[j]][[i]])
		}
		highOrderNbList$total[[i]] <- sort(indx)
	}
	
	return(highOrderNbList)
}
