
split.predictions.createDistMat <- function(ssn, predpts = NULL, o.write = FALSE)
	{

	if(missing(predpts) || is.null(predpts) || length(predpts) == 0)
	{
		stop("A named collection of prediction points must be specified via the predpts option")
	}
	#Check to see whether distance folder exists...
	if (!file.exists(file.path(ssn@path, "distance")))
	{
		dir.create(file.path(ssn@path, "distance"))
	}
	#And then whether relevant predpts folders exist
	for(pred in predpts)
	{
		if(!file.exists(file.path(ssn@path, "distance", pred)))
		{
			dir.create(file.path(ssn@path, "distance", pred))
		}
    }
	#cast netID to a factor for every set of predpoints
	for(i in 1:length(ssn@predpoints@SSNPoints))
	{
		ssn@predpoints@SSNPoints[[i]]@point.data$netID<-
                    as.factor(ssn@predpoints@SSNPoints[[i]]@point.data$netID)
	}
	#and also for the observed points
	ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID<-
            as.factor(ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID)

	#Initialise binaryID.db
	if (file.exists(file.path(ssn@path,"binaryID.db")) == FALSE)
	{
		stop("binaryID.db is missing from ssn object")
	}

	driver <- RSQLite::SQLite()
	connect.name <- file.path(ssn@path,"binaryID.db")
	connect <- dbConnect(SQLite(), connect.name)
	#close sqlite connection on function exit
	on.exit({
          dbDisconnect(connect)
	  ##sqliteCloseConnection(connect)
	  ##sqliteCloseDriver(driver)
	}, add=TRUE)

	net.count <- length(levels(ssn@network.line.coords$NetworkID))
	warned.overwrite <- FALSE
	for(pred.num in 1:length(ssn@predpoints@SSNPoints))
	{
		#If we're not interested in this set of prediction points, skip it
		if(!(ssn@predpoints@ID[pred.num] %in% predpts)) next
		predpt.name <- ssn@predpoints@ID[pred.num]
		for (i in 1:net.count)
		{
			#get the number (rather than index) of this network
			net.num <- levels(ssn@network.line.coords$NetworkID)[i]

			#find out which prediction points of this set are on the specified network
			ind.preds <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$NetworkID == as.numeric(net.num)
			pred.site.no <- sum(ind.preds)

			#if there are none continue
			if(pred.site.no == 0) next

			#figure out how many observed sites there are in the network
			ind.obs <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == as.numeric(net.num)
			site.no <- nrow(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])

			# get sorted pids to use as dim names
			obs.pids <- sort(as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])))
			pred.pids<- sort(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])))

			net.name <- paste("net", net.num, sep = "")
			workspace.name.a <- paste("dist.net", net.num, ".a.RData", sep = "")
			workspace.name.b <- paste("dist.net", net.num, ".b.RData", sep = "")

			bin.table <- dbReadTable(connect, net.name)
			exists <- c(file.exists(file.path(ssn@path, "distance", predpt.name, workspace.name.a)),
                                    file.exists(file.path(ssn@path, "distance", predpt.name, workspace.name.b)))
			if(!o.write)
			{
				if(all(exists))
				{
					if(!warned.overwrite)
					{
						warned.overwrite <- TRUE
						##warning("One or more set of observed - prediction distance matrices already existed and was not overwritten")
					}
					next
				}
			}
			current_distance_matrix_a <- matrix(NA, nrow = site.no, ncol = pred.site.no, dimnames=list(obs.pids, pred.pids))
			current_distance_matrix_b <- matrix(NA, nrow = pred.site.no, ncol = site.no, dimnames=list(pred.pids, obs.pids))

			if (site.no !=0)
			{
				#Start of incomprehensible chunk lifted from createDistMat code.
				locID.obi <- attributes(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])$locID
				ob.i <- as.data.frame(cbind(as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])),
				as.numeric(levels(ssn@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[ind.obs]))[ssn@obspoints@SSNPoints[[1]]@network.point.coords$SegmentID[ind.obs]],
				locID.obi[ind.obs]))

				colnames(ob.i)<- c("pid","rid","locID")
				ob.i$locID <- as.factor(ob.i$locID)
				ob.i$binaryID <- bin.table$binaryID[match(ob.i$rid, bin.table$rid)]
				ob.i <-ob.i[order(ob.i[,"pid"]),]
				rownames(ob.i)<- ob.i$pid

				ob.i_by_locID <- ob.i[order(ob.i[,"locID"]),]
				ob.i_by_locID$pid <- as.numeric(ob.i_by_locID$pid)
				ob.i_by_locID$locID <- as.numeric(ob.i_by_locID$locID)
				ob.j_reordering <- order(ob.i_by_locID$pid)

				locID.old <- -1
				ind.dup <- !duplicated(ob.i_by_locID$locID)

				for (j in 1:nrow(ob.i))
				{
					pid.i <- ob.i[j,"pid"]
					locID.i <- ob.i[j, "locID"]

					if (locID.i != locID.old)
					{
						upDist.i <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[paste(pid.i),"DistanceUpstream"]

						pred.tmp <- as.data.frame(cbind(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])),
                                                     as.numeric(levels(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[ind.preds]))[ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[ind.preds]]))
						colnames(pred.tmp)<- c("pid","rid")

						pred.tmp$binaryID <- bin.table$binaryID[match(pred.tmp$rid, bin.table$rid)]
						pred.tmp <-pred.tmp[order(pred.tmp[,"pid"]),]
						rownames(pred.tmp) <- pred.tmp$pid

						junk <- get.rid.fc(pred.tmp[,"binaryID"], ob.i$binaryID[j])
						ob.j <- data.frame(pred.tmp["pid"], junk, stringsAsFactors = FALSE)

						ob.j$pid <- as.numeric(ob.j$pid)
						ob.j$fc <- as.logical(ob.j$fc)

						ob.j$junc.rid <- bin.table$rid[match(ob.j$binaryID, bin.table$binaryID)]

						ob.j$juncDist <- ssn@network.line.coords$DistanceUpstream[match(ob.j$junc.rid, ssn@network.line.coords$SegmentID)]

						ob.j$upDist.j <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$DistanceUpstream[match(ob.j$pid,
                                                                 as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords)))]

						ob.j <-ob.j[order(ob.j[,"pid"]),]

						ind.fc<-ob.j$fc==1
						dist.a <- ifelse(ind.fc, ob.j$upDist.j-upDist.i, ob.j$upDist.j - ob.j$juncDist)
						current_distance_matrix_a[paste(pid.i), ] <- ifelse(dist.a<0, 0, dist.a)

						dist.b <- ifelse(ind.fc, upDist.i - ob.j$upDist.j, upDist.i - ob.j$juncDist)
						current_distance_matrix_b[, paste(pid.i)] <- ifelse(dist.b<0, 0, dist.b)
					}
					else
					{
						# add column to pred sites
						current_distance_matrix_a[paste(pid.i),]<- current_distance_matrix_a[paste(pid.old),]
						current_distance_matrix_b[,paste(pid.i)]<- current_distance_matrix_b[,paste(pid.old)]
					}
				}
				pid.old <- pid.i
				locID.old <- locID.i
			}
			#save obs-pred and pred-obs distance matrices
			file_handle = file(file.path(ssn@path, "distance", predpt.name, workspace.name.a), open="wb")
			serialize(current_distance_matrix_a, file_handle, ascii=FALSE)
			close(file_handle)

			file_handle = file(file.path(ssn@path, "distance", predpt.name, workspace.name.b), open="wb")
			serialize(current_distance_matrix_b, file_handle, ascii=FALSE)
			close(file_handle)
		}
	}
}

splitPredictions <- function(ssn, predpointsID, chunksof, by, subset, new.id)
{
	if(missing(ssn) || missing(predpointsID)) stop("Arguments ssn and predpointsID must be specified")

	if(sum(c(!missing(chunksof), !missing(by), !missing(subset))) != 1) stop("Only one of 'chunksof', 'by' and 'subset' can be input to splitPredictions")
	if(!missing(new.id) && missing(subset)) stop("Input 'new.id' should only be specified if input 'subset' is used")
	if(!missing(subset) && missing(new.id)) stop("Input 'new.id' must be specified if input 'subset' is used")

	pred.len <- length(ssn@predpoints@ID)
	found.pred <- FALSE
	for(i in 1:pred.len)
	{
		if(ssn@predpoints@ID[[i]] == predpointsID)
		{
			point.coords.size <- nrow(ssn@predpoints@SSNPoints[[i]]@point.coords)
			point.data.size <- nrow(ssn@predpoints@SSNPoints[[i]]@point.data)
			network.point.coords.size <- nrow(ssn@predpoints@SSNPoints[[i]]@network.point.coords)
			if(point.coords.size != point.data.size || point.data.size != network.point.coords.size) stop("Size mismatches")

                        pid.vec <- ssn@obspoints@SSNPoints[[1]]@point.data$pid
                        for(z in 1:pred.len) {
                            pid.vec <- append(pid.vec, ssn@predpoints@SSNPoints[[z]]@point.data$pid)
                        }
                        if(sum(duplicated(pid.vec))>0) stop("Duplicated pid values in ssn")
                        max.pid <- max(pid.vec)

			old.wd <- getwd()
			on.exit(setwd(old.wd), add=TRUE)
			setwd(ssn@path)

			if(!missing(chunksof))
			{
				npoints <- point.data.size
				nchunks <- round((npoints / chunksof) + 0.5)
				if(nchunks > 1000) stop(paste("Specified value of chunksof would result in prediction points set ",
                                        predpointsID, " being split into more than 1000 prediction sets", sep=""))
				chunks <- seq.int(1, point.data.size, chunksof)
				if(tail(chunks, 1) != point.data.size) chunks <- c(chunks, point.data.size+1)
				for(j in 1:(length(chunks)-1))
				{
					new.id <- paste(predpointsID, "-", j, sep="")
					subsetted.coords <- ssn@predpoints@SSNPoints[[i]]@point.coords[chunks[j]:(chunks[j+1]-1), ,drop=FALSE]
					subsetted.data <- ssn@predpoints@SSNPoints[[i]]@point.data[chunks[j]:(chunks[j+1]-1), ,drop=FALSE]
					subsetted.network.point.coords <- ssn@predpoints@SSNPoints[[i]]@network.point.coords[chunks[j]:(chunks[j+1]-1), ,drop=FALSE]

                                        ## Set subsetted.coords pid
                                        ##sub.size <- nrow(subsetted.coords)
                                        tmp.pid <- as.numeric(row.names(subsetted.coords))
                                        tmp.pid <- tmp.pid + max.pid
                                        if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                        rownames(subsetted.coords) <-  format(tmp.pid, scientific = FALSE)
                                        rm(tmp.pid)

                                        ## set subsetted.data pid
                                        if(sum(subsetted.data$pid != rownames(subsetted.data))>0) stop("rownames do not match pid values in point.data")
                                        tmp.pid <- subsetted.data$pid + max.pid
                                        if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                        rownames(subsetted.data) <- tmp.pid
                                        subsetted.data$pid <- tmp.pid
                                        rm(tmp.pid)

                                        ## set subsetted.network.point.coords pid
                                        tmp.pid <- as.numeric(row.names(subsetted.network.point.coords))
                                        tmp.pid <- tmp.pid + max.pid
                                        if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                        rownames(subsetted.network.point.coords) <-  format(tmp.pid, scientific = FALSE)
                                        attributes(subsetted.network.point.coords)$locID <- as.numeric(levels(subsetted.data$locID))[subsetted.data$locID]
                                        pid.vec <- append(pid.vec, tmp.pid)
                                        rm(tmp.pid)

					##write to file
					subsetted.spatialStructure <- sp::SpatialPointsDataFrame(coords = subsetted.coords,
                                                  data = subsetted.data, proj4string = ssn@predpoints@SSNPoints[[i]]@proj4string)
					maptools::writeSpatialShape(subsetted.spatialStructure, new.id)

					#and alter the existing object
					new.index <- length(ssn@predpoints@ID) + 1
					ssn@predpoints@ID[[new.index]] <- new.id
					ssn@predpoints@SSNPoints[[new.index]] <- new("SSNPoint", point.coords = subsetted.coords,
                                                  point.data = subsetted.data, network.point.coords = subsetted.network.point.coords,
                                                  points.bbox = ssn@predpoints@SSNPoints[[i]]@points.bbox,
                                                  proj4string = ssn@predpoints@SSNPoints[[i]]@proj4string)
				}
				new.predids <- paste(predpointsID,"-", 1:(length(chunks)-1), sep="")
			}
			else if(!missing(subset))
			{
				e <- substitute(subset)
				values <- eval(e, ssn@predpoints@SSNPoints[[i]]@point.data, parent.frame())
				if(!is.logical(values)) stop("Input 'subset' must evaluate to a logical value")

				subsetted.coords <- ssn@predpoints@SSNPoints[[i]]@point.coords[values, ,drop=FALSE]
				subsetted.data <- ssn@predpoints@SSNPoints[[i]]@point.data[values, ,drop=FALSE]
				subsetted.network.point.coords <- ssn@predpoints@SSNPoints[[i]]@network.point.coords[values, ,drop=FALSE]

                                ## Set subsetted.coords pid
                                tmp.pid <- as.numeric(row.names(subsetted.coords))
                                tmp.pid <- tmp.pid + max.pid
                                if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                rownames(subsetted.coords) <-  format(tmp.pid, scientific = FALSE)
                                rm(tmp.pid)

                                ## set subsetted.data pid
                                if(sum(subsetted.data$pid != rownames(subsetted.data))>0) stop("rownames do not match pid values in point.data")
                                tmp.pid <- subsetted.data$pid + max.pid
                                if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                rownames(subsetted.data) <- tmp.pid
                                subsetted.data$pid <- tmp.pid
                                rm(tmp.pid)

                                ## set subsetted.network.point.coords pid
                                tmp.pid <- as.numeric(row.names(subsetted.network.point.coords))
                                tmp.pid <- tmp.pid + max.pid
                                if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                rownames(subsetted.network.point.coords) <-  format(tmp.pid, scientific = FALSE)
                                attributes(subsetted.network.point.coords)$locID <- as.numeric(levels(subsetted.data$locID))[subsetted.data$locID]
                                pid.vec <- append(pid.vec, tmp.pid)
                                rm(tmp.pid)

				subsetted.spatialStructure <- sp::SpatialPointsDataFrame(coords = subsetted.coords,
                                          data = subsetted.data, proj4string = ssn@predpoints@SSNPoints[[i]]@proj4string)
				maptools::writeSpatialShape(subsetted.spatialStructure, new.id)

				#and alter the existing object
				new.index <- length(ssn@predpoints@ID) + 1
				ssn@predpoints@ID[[new.index]] <- new.id

				ssn@predpoints@SSNPoints[[new.index]] <- new("SSNPoint", point.coords = subsetted.coords,
                                      point.data = subsetted.data, network.point.coords = subsetted.network.point.coords,
                                      points.bbox = ssn@predpoints@SSNPoints[[i]]@points.bbox,
                                      proj4string = ssn@predpoints@SSNPoints[[i]]@proj4string)
				new.predids <- new.id
			}
			else
			{
				if(!(by %in% colnames(ssn@predpoints@SSNPoints[[i]]@point.data))) stop(paste("Could not find column named ",by, " in point.coords entry of specified prediction points set", sep=""))
				current.class <- class(ssn@predpoints@SSNPoints[[i]]@point.data[[by]])
				if(current.class != "factor") warning(paste("Column named ", by, " had class ", current.class, " instead of class factor. Casting to factor....", sep=""))
				ssn@predpoints@SSNPoints[[i]]@point.data[[by]] <- as.factor(ssn@predpoints@SSNPoints[[i]]@point.data[[by]])
				levels <- levels(ssn@predpoints@SSNPoints[[i]]@point.data[[by]])
				new.predids <- c()
				for(level in levels)
				{
					new.id <- paste(predpointsID, "-", by, "-", level, sep="")
					new.predids <- c(new.predids, new.id)

					relevant <- ssn@predpoints@SSNPoints[[i]]@point.data[[by]] == level
					subsetted.data <- ssn@predpoints@SSNPoints[[i]]@point.data[relevant,, drop=FALSE]
					subsetted.coords <- ssn@predpoints@SSNPoints[[i]]@point.coords[relevant,, drop=FALSE]
					subsetted.network.point.coords <- ssn@predpoints@SSNPoints[[i]]@network.point.coords[relevant,,drop=FALSE]
					subsetted.data[[by]] <- rep(level, length(subsetted.data[[by]]))

                                        ## Set subsetted.coords pid
                                        tmp.pid <- as.numeric(row.names(subsetted.coords))
                                        tmp.pid <- tmp.pid + max.pid
                                        if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                        rownames(subsetted.coords) <-  format(tmp.pid, scientific = FALSE)
                                        rm(tmp.pid)

                                        ## set subsetted.data pid
                                        if(sum(subsetted.data$pid != rownames(subsetted.data))>0) stop("rownames do not match pid values in point.data")
                                        tmp.pid <- subsetted.data$pid + max.pid
                                        if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                        rownames(subsetted.data) <- tmp.pid
                                        subsetted.data$pid <- tmp.pid
                                        rm(tmp.pid)

                                        ## set subsetted.network.point.coords pid
                                        tmp.pid <- as.numeric(row.names(subsetted.network.point.coords))
                                        tmp.pid <- tmp.pid + max.pid
                                        if(sum(tmp.pid %in% pid.vec)>0) stop("Duplicate pids exist")
                                        rownames(subsetted.network.point.coords) <- format(tmp.pid, scientific = FALSE)
                                        attributes(subsetted.network.point.coords)$locID <- as.numeric(levels(subsetted.data$locID))[subsetted.data$locID]
                                        pid.vec <- append(pid.vec, tmp.pid)
                                        rm(tmp.pid)

					#write to file
					subsetted.spatialStructure <- sp::SpatialPointsDataFrame(coords = subsetted.coords,
                                            data = subsetted.data, proj4string = ssn@predpoints@SSNPoints[[i]]@proj4string)
					maptools::writeSpatialShape(subsetted.spatialStructure, new.id)

					#and alter the existing object
					new.index <- length(ssn@predpoints@ID) + 1
					ssn@predpoints@ID[[new.index]] <- new.id

					ssn@predpoints@SSNPoints[[new.index]] <- new("SSNPoint", point.coords = subsetted.coords,
                                            point.data = subsetted.data, network.point.coords = subsetted.network.point.coords,
                                            points.bbox = ssn@predpoints@SSNPoints[[i]]@points.bbox,
                                            proj4string = ssn@predpoints@SSNPoints[[i]]@proj4string)
				}
			}
			found.pred <- TRUE
		}

	}
	if(!found.pred) stop(paste("Given object did not contain prediction points ID ", predpointsID, sep=""))

	split.predictions.createDistMat(ssn, new.predids, o.write=FALSE)
	return(ssn)
}
