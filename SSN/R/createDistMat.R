createDistMat <-
function(ssn, predpts = NULL, o.write = FALSE, amongpreds = FALSE) {

  #start.time <- Sys.time()
  if(amongpreds && (missing(predpts) || is.null(predpts)))
  {
	stop("A named collection of prediction points must be specified via the predpts option when amongpreds is TRUE")
  }
  ##Check to see whether distance folder exists...
  if (!file.exists(file.path(ssn@path, "distance"))) {
    dir.create(file.path(ssn@path, "distance"))
  }
  ##And then whether an observation folder exists
  if (!file.exists(file.path(ssn@path, "distance", "obs"))) {
    dir.create(file.path(ssn@path, "distance", "obs"))
  }

  ## And then whether prediction folder exists
  if (!is.null(predpts)) {
    if(!file.exists(file.path(ssn@path, "distance", predpts))) {
      dir.create(file.path(ssn@path, "distance", predpts))
    }
    count <- 0
    if(length(ssn@predpoints@ID) > 0) {
	for (m in 1:length(ssn@predpoints@ID)) {
	    if (ssn@predpoints@ID[m] == predpts) {
	         pred.num <- m
			count <- count + 1}
		}
	}

    if (count==0) {
      stop(predpts, " does not exist in SSN")
    }
    if (count > 1) {
      stop("SSN contains more than one copy of ", predpts)
    }
	ssn@predpoints@SSNPoints[[pred.num]]@point.data$netID<- as.factor(ssn@predpoints@SSNPoints[[pred.num]]@point.data$netID)
  }
  ##########
  if (is.null(predpts)) {
      pred.num <- 0
  }

  ##Initialise binaryID.db
  if (file.exists(file.path(ssn@path,"binaryID.db")) == FALSE)
    stop("binaryID.db is missing from ssn object")

  driver <- RSQLite::SQLite()
  connect.name <- file.path(ssn@path,"binaryID.db")

  connect <- dbConnect(SQLite(), connect.name)

  ## close sqlite connection upon function exit
  on.exit({
      dbDisconnect(connect)
      ##sqliteCloseConnection(connect)
      ##sqliteCloseDriver(driver)
  })

  if (file.exists(file.path(ssn@path, "binaryID.db")) == TRUE) {

    ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID<- as.factor(ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID)
    net.count <- length(levels(ssn@network.line.coords$NetworkID))
    warned.overwrite <- FALSE
    for (i in 1:net.count) {
      net.num <- levels(ssn@network.line.coords$NetworkID)[i]
      ind.obs <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == as.numeric(net.num)

      ## figure out how many observed and prediction sites there are in the network
      site.no <- nrow(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])


      ####
      if (pred.num > 0) {
          ind.preds <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$NetworkID == as.numeric(net.num)
          pred.site.no <- nrow(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])
      } else { pred.site.no <-0 }

      if (site.no > 0) {

      ## get sorted pids to use as dim names
      ##obs.pids <- sort(as.numeric(ssn@obspoints@SSNPoints[[1]]@point.data[ind.obs,"pid"]))
	  obs.pids<- sort(as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.obs,])))


####################################################
## If predpts !is.null
####################################################

	  if(!is.null(predpts)){

		##ind.preds <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$NetworkID == as.numeric(net.num)
		##pred.site.no <- nrow(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])
		pred.pids<- sort(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])))



####################################################
## if pred.site.no > 0
####################################################

		if(pred.site.no > 0){

		    ## create o x p distance matrix full of NA
                    current_distance_matrix_a <- matrix(NA, nrow = site.no, ncol = pred.site.no,
                    dimnames=list(obs.pids, pred.pids))

                    ## create p x o distance matrix full of NA
                    current_distance_matrix_b <- matrix(NA, nrow = pred.site.no, ncol = site.no,
                    dimnames=list(pred.pids, obs.pids))
		}

	  }

          net.name <- paste("net", net.num, sep = "")
          workspace.name.a <- paste("dist.net", net.num, ".a.RData", sep = "")
	  workspace.name.b <- paste("dist.net", net.num, ".b.RData", sep = "")

          bin.table <- dbReadTable(connect, net.name)

  	  ##Create n x n distance matrix full of NA
	  ##If distance matrix table already exists within the database, exit function
	  workspace.name <- paste("dist.net", net.num, ".RData", sep = "")

	  if(!o.write) {
		exists <- file.exists(file.path(ssn@path, "distance", "obs", workspace.name))
		if (!missing(predpts) && !is.null(predpts))
		{
			exists <- c(exists, file.exists(file.path(ssn@path, "distance", predpts, workspace.name.a)), file.exists(file.path(ssn@path, "distance", predpts, workspace.name.b)))
		}
		##If they're all already there, warn but continue
		if(all(exists))
		{
			if(!warned.overwrite)
			{
				warned.overwrite <- TRUE
				cat("Distance matrices already existed while o.write was set to FALSE. Not overwriting existing matrices\n")
			}
			next
		}
		##if some exist but some don't that's an error
		else if(any(exists) && any(!exists))
		{
			stop("o.write was set to FALSE and some (but not all) distance matrices already existed")
		}
	  }

          ## Create obs distance matrix
	  current_distance_matrix <- matrix(NA, nrow = site.no, ncol = site.no,dimnames = list(obs.pids, obs.pids))
	  diag(current_distance_matrix)<- 0

          rownames(current_distance_matrix) <- obs.pids
	  colnames(current_distance_matrix) <- obs.pids

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

##########Calculate distances between observed sites
          for (j in 1:nrow(ob.i)) {

            pid.i <- ob.i[j,"pid"]
            locID.i <- ob.i[j, "locID"]

            if (locID.i != locID.old) {
	      upDist.i <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[paste(pid.i),"DistanceUpstream"]
              junk <- get.rid.fc(ob.i_by_locID[ind.dup,"binaryID"], ob.i$binaryID[j])

              ob.j.r <- data.frame(ob.i_by_locID[ind.dup, c("pid", "locID")], junk, stringsAsFactors = FALSE)

              ob.j.r$fc <- as.logical(ob.j.r$fc)
              rownames(ob.j.r)<- ob.j.r$pid

              ob.j.r$junc.rid <- bin.table$rid[match(ob.j.r$binaryID, bin.table$binaryID)]

              reps <- as.numeric(ob.i_by_locID$locID)
              ob.j <- ob.j.r[reps,]

              rownames(ob.j) <- paste(rownames(ob.i), ".fc", sep = "")

              ob.j$pid <- ob.i_by_locID$pid
	      ob.j$juncDist <- ssn@network.line.coords$DistanceUpstream[match(ob.j$junc.rid, ssn@network.line.coords$SegmentID)]

	      ob.j$upDist.j <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream[
              match(ob.j$pid, as.numeric(rownames(ssn@obspoints@SSNPoints[[1]]@network.point.coords)))]

              ob.j <-ob.j[ob.j_reordering,]


              ##obs fills in by column because it's working between obs to obs.
              ind.fc<-ob.j$fc==1

              dist.obs <- ifelse(ind.fc, upDist.i - ob.j$upDist.j, upDist.i - ob.j$juncDist)
              current_distance_matrix[,paste(pid.i)] <- ifelse(dist.obs<0, 0, dist.obs)


    #####Pred sites--------------------------------------------------------------
          } else {
              current_distance_matrix[,paste(pid.i)]<- current_distance_matrix[,paste(pid.old)]
          }

          if (locID.i != locID.old) {

            if (!is.null(predpts) && pred.site.no > 0) {

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
	      ob.j$upDist.j <- ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$DistanceUpstream[
              match(ob.j$pid, as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords)))]

              ob.j <-ob.j[order(ob.j[,"pid"]),]

              ind.fc<-ob.j$fc==1
              dist.a <- ifelse(ind.fc, ob.j$upDist.j-upDist.i, ob.j$upDist.j - ob.j$juncDist)
              current_distance_matrix_a[paste(pid.i), ] <- ifelse(dist.a<0, 0, dist.a)

              dist.b <- ifelse(ind.fc, upDist.i - ob.j$upDist.j, upDist.i - ob.j$juncDist)
              current_distance_matrix_b[, paste(pid.i)] <- ifelse(dist.b<0, 0, dist.b)
            }
          } else {
            ## add column to pred sites
            if (!is.null(predpts)) {
              current_distance_matrix_a[paste(pid.i),]<- current_distance_matrix_a[paste(pid.old),]
              current_distance_matrix_b[,paste(pid.i)]<- current_distance_matrix_b[,paste(pid.old)]}
          }

          pid.old <- pid.i
          locID.old <- locID.i
        }

      ##save distance matrix
      file_handle = file(file.path(ssn@path, "distance", "obs", workspace.name), open="wb")
      serialize(current_distance_matrix, file_handle, ascii=FALSE)
      close(file_handle)


      if(pred.site.no > 0) {
      	##save obs-pred and pred-obs distance matrices
      	file_handle = file(file.path(ssn@path, "distance", predpts, workspace.name.a), open="wb")
      	serialize(current_distance_matrix_a, file_handle, ascii=FALSE)
      	close(file_handle)

      	file_handle = file(file.path(ssn@path, "distance", predpts, workspace.name.b), open="wb")
      	serialize(current_distance_matrix_b, file_handle, ascii=FALSE)
      	close(file_handle)

          ##if(amongpreds && !is.null(predpts)) {
          if(amongpreds) {
             ##Create distance matrix full of NAs. If the file already exists and overwrite == FALSE, stop with an error
             workspace.name <- paste("dist.net", net.num, ".RData", sep = "")

             among_distance_matrix <- matrix(NA, nrow = pred.site.no, ncol = pred.site.no)
             diag(among_distance_matrix)<- 0

             rownames(among_distance_matrix) <- pred.pids
             colnames(among_distance_matrix) <- pred.pids


	    locID.pred.data <- attributes(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords)$locID
            pred.data <- as.data.frame(cbind(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])),
				as.numeric(levels(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[ind.preds]))[ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[ind.preds]],
                                locID.pred.data[ind.preds], ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$DistanceUpstream[ind.preds]))

            colnames(pred.data)<- c("pid","rid", "locID", "upDist")
            pred.data <- pred.data[order(pred.data$pid),]
            pred.data$binaryID <- bin.table$binaryID[match(pred.data$rid, bin.table$rid)]
	    pred.data <- pred.data[order(pred.data[,"pid"]),]
            rownames(pred.data) <- pred.data$pid

            ##locID values can be repeated, in which case they have the same distance data.
	    locID.old <- -1
	    for(b in 1:pred.site.no){
	       locID.b <- pred.data[b, "locID"]
	       upDist.b <- pred.data[b, "upDist"]
	       pid.b <- pred.data[b, "pid"]

               if(locID.b != locID.old) {
                   junk <- get.rid.fc(pred.data[,"binaryID"], pred.data$binaryID[b])
                   truncated.binaryIDs <- data.frame(pid=pred.data[,"pid"], junk, stringsAsFactors = FALSE)
                   truncated.binaryIDs$fc <- as.logical(truncated.binaryIDs$fc)
                   truncated.binaryIDs$junc.rid <- bin.table$rid[match(truncated.binaryIDs$binaryID, bin.table$binaryID)]

                   truncated.binaryIDs$juncDist <- ssn@data$upDist[match(truncated.binaryIDs$junc.rid, ssn@data$rid)]
                   truncated.binaryIDs$upDist.j <- pred.data$upDist[match(truncated.binaryIDs$pid, pred.data$pid)]
                   ind.fc<-truncated.binaryIDs$fc==1
                   dist.preds <- ifelse(ind.fc, upDist.b - truncated.binaryIDs$upDist.j, upDist.b - truncated.binaryIDs$juncDist)
                   among_distance_matrix[,paste(pid.b)] <- ifelse(dist.preds<0, 0, dist.preds)
	       }
	    }

             file_handle = file(file.path(ssn@path, "distance", predpts, workspace.name), open="wb")
	    serialize(among_distance_matrix, file_handle, ascii=FALSE)
	    close(file_handle)
      }}
      }

      if (site.no ==0 & pred.site.no > 0) {
          if(amongpreds) {
             ##Create distance matrix full of NAs. If the file already exists and overwrite == FALSE, stop with an error
             workspace.name <- paste("dist.net", net.num, ".RData", sep = "")
             pred.pids<- sort(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])))
             net.name <- paste("net", net.num, sep = "")
             bin.table <- dbReadTable(connect, net.name)

             among_distance_matrix <- matrix(NA, nrow = pred.site.no, ncol = pred.site.no)
             diag(among_distance_matrix)<- 0

## Add pred.pids into this function
####################################

             rownames(among_distance_matrix) <- pred.pids
             colnames(among_distance_matrix) <- pred.pids


	    locID.pred.data <- attributes(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords)$locID
            pred.data <- as.data.frame(cbind(as.numeric(rownames(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords[ind.preds,])),
				as.numeric(levels(ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[ind.preds]))[ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$SegmentID[ind.preds]],
                                locID.pred.data[ind.preds], ssn@predpoints@SSNPoints[[pred.num]]@network.point.coords$DistanceUpstream[ind.preds]))

            colnames(pred.data)<- c("pid","rid", "locID", "upDist")
            pred.data <- pred.data[order(pred.data$pid),]

## Need bin.table

            pred.data$binaryID <- bin.table$binaryID[match(pred.data$rid, bin.table$rid)]
	    pred.data <- pred.data[order(pred.data[,"pid"]),]
            rownames(pred.data) <- pred.data$pid

            ##locID values can be repeated, in which case they have the same distance data.
	    locID.old <- -1
	    for(b in 1:pred.site.no){
	       locID.b <- pred.data[b, "locID"]
	       upDist.b <- pred.data[b, "upDist"]
	       pid.b <- pred.data[b, "pid"]

               if(locID.b != locID.old) {
                   junk <- get.rid.fc(pred.data[,"binaryID"], pred.data$binaryID[b])
                   truncated.binaryIDs <- data.frame(pid=pred.data[,"pid"], junk, stringsAsFactors = FALSE)
                   truncated.binaryIDs$fc <- as.logical(truncated.binaryIDs$fc)
                   truncated.binaryIDs$junc.rid <- bin.table$rid[match(truncated.binaryIDs$binaryID, bin.table$binaryID)]

                   truncated.binaryIDs$juncDist <- ssn@data$upDist[match(truncated.binaryIDs$junc.rid, ssn@data$rid)]
                   truncated.binaryIDs$upDist.j <- pred.data$upDist[match(truncated.binaryIDs$pid, pred.data$pid)]
                   ind.fc<-truncated.binaryIDs$fc==1
                   dist.preds <- ifelse(ind.fc, upDist.b - truncated.binaryIDs$upDist.j, upDist.b - truncated.binaryIDs$juncDist)
                   among_distance_matrix[,paste(pid.b)] <- ifelse(dist.preds<0, 0, dist.preds)
	       }
	    }

            file_handle = file(file.path(ssn@path, "distance", predpts, workspace.name), open="wb")
	    serialize(among_distance_matrix, file_handle, ascii=FALSE)
	    close(file_handle)
      }}
}}}

