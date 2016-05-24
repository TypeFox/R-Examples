construct.s <- function(spikes, ids, time.interval, beg, end, corr.breaks, 
                        layout, filename) {
  spikes.range <- range(unlist(spikes))
  if (!is.null(end)) {
    spikes <- lapply(spikes, .jay.filter.for.max, max = end)
  } else {
    end <- spikes.range[2]
  }
  if (!is.null(beg)) {
    spikes <- lapply(spikes, .jay.filter.for.min, min = beg)
  } else {
    beg <- spikes.range[1]
  }
  
  ## Remove any channels that have zero spikes (which can happen if
  ## the beg, end range is too narrow, or if a datafile is empty,
  ## which happens sometime for the Feller data.
  spikes <- .remove.empty.channels(spikes)
  if (!is.null(ids)) {
    spikes <- .filter.channel.names(spikes, ids)
  }
  
  ## No more spike trains should be removed, so we can refilter the layout to ensure
  ## that only spikes that are left in the array are kept in the layout.
  channels <- names(spikes)
  keep <- match(channels, rownames(layout$pos))
  layout$pos <- layout$pos[keep, ]
  rec.time <- c(beg, end)
  nspikes <- sapply(spikes, length)
  
  ###################### patch
  if (length(nspikes) ==0) {
    meanfiringrate <- nspikes #empty list
  } else {
    meanfiringrate <- nspikes/(end - beg)
  }
  ######################
  
  rates <- .make.spikes.to.frate(spikes, time.interval = time.interval, 
                                beg = beg, end = end)
  unit.offsets <- NULL
  .check.spikes.monotonic(spikes)
  res <- list(channels = names(spikes), spikes = spikes, nspikes = nspikes, 
              NCells = length(spikes), meanfiringrate = meanfiringrate, 
              file = filename, layout = layout, rates = rates, unit.offsets = unit.offsets, 
              rec.time = rec.time)
  class(res) <- "mm.s"
  if (length(corr.breaks) == 1) {
    res$corr = NULL
  }    else {
    res$corr = .corr.index(res, corr.breaks)
  }
  res
}



calculate.isis <- function(s) {
	
	s$isis <- lapply(s$spikes,diff)
	s$mean.isis <- lapply(s$isis,mean)
	s$sd.isis <- lapply(s$isis,sd)
	return(s)
}

.plot.isis.by.plate <- function(s) {
	isis.all <- unlist(s$isis)
	hist(isis.all,main = "Histogram of ISIs by Plate", xlab = "ISI length")
	hist(log10(isis.all),main = "Histogram of log(ISIs) by Plate", xlab = "log10(ISI length)")
}
.spike.summary.by.electrode <- function(s) {
	s <- calculate.isis(s)
	electrodes <- .get.all.electrodes(s)
	sum <- matrix(data = NA, nrow = length(electrodes), ncol = 4)
	colnames(sum) <- c("nspikes","meanfiringrate","meanisis","sdisis")
	rownames(sum) <- electrodes

	df <- cbind(s$nspikes,s$meanfiringrate,s$mean.isis,s$sd.isis)
	active.electrodes<-rownames(df)
	for (i in active.electrodes) {
		sum[i,] <- unlist(df[i,])
	}
	sum
}

.spike.summary.by.well <- function(s) {
  plate <- plateinfo(s$layout$array)
  wells <- sort(plate$wells)
  s$isis <- lapply(s$spikes,diff)
  # Sahar 10282014 - add genotype in 1st column to output and generic column position (startPos) to add more columns
  startPos=1
  sum <- matrix(data = NA, nrow = length(wells), ncol = startPos+8)
  colnames(sum)<- c("treatment","nAE","nspikes_by_well","meanfiringrate_by_well","meanfiringrate_by_all_ectctordes","meanfiringrate_by_active_electordes","sdfiringrate_by_active_electordes","meanisis","sdisis")
  rownames(sum) <- wells
  nelectrodes <- plate$n.elec.r *plate$n.elec.c
  if (!is.null(s$goodwells)) {
    for (j in 1:length(s$goodwells)){
      icurrentwell<-(s$goodwells[j]==s$cw)
      incurrentwell<-which((s$goodwells[j]==s$cw))
      #sahar - add genotype from treatment column that is per well (if well has data), not per electrode
      #treatment="NA"
      # diana change because trt was showing as na when there was a trt, just no activity
      # check that current well and treatment well match

      treatment=s$treatment[s$goodwells[j]]
      
      if (length(incurrentwell)>0){
        well <- strsplit(s$channels[incurrentwell], "_")[[1]][1]
        treatment=s$treatment[well][[1]] }
      sum[s$goodwells[j],startPos] <- treatment
      sum[s$goodwells[j],startPos+1] <- length(incurrentwell)
      sum[s$goodwells[j],startPos+2] <- sum(s$nspikes[icurrentwell])
      sum[s$goodwells[j],startPos+3] <- sum(s$meanfiringrate[icurrentwell])
      sum[s$goodwells[j],startPos+4] <- sum(s$meanfiringrate[icurrentwell])/ nelectrodes #Sahar - was division of string by number - ask Quanli!
      
      sum[s$goodwells[j],startPos+5] <- mean(s$meanfiringrate[icurrentwell])
      sum[s$goodwells[j],startPos+6] <- sd(s$meanfiringrate[icurrentwell])
      
      isis.all <- unlist(s$isis[icurrentwell])
      sum[s$goodwells[j],startPos+7] <- mean(isis.all)
      sum[s$goodwells[j],startPos+8] <- sd(isis.all)
    }
  }
  sum
}



.get.div<-function (s) {
  div <- NA
  t1 <- strsplit(s$file, split="_", fixed = TRUE)
  for (i in t1[[1]]) {
    if (nchar(i) > 2 && substr(i, 1, 3) == "DIV") {
      if (nchar(i)>5){
        i=unlist(strsplit(i, split=".", fixed=T))[1]
      }
      div <- as.numeric(substr(i, 4, nchar(i)))
    }
  }
  div
}

IGM.mean.firingrate.by.well <- function(s) {
	df1 <- aggregate(s$meanfiringrate, by = list(s$cw),FUN = mean,na.rm = T)
	df2 <- aggregate(s$meanfiringrate, by = list(s$cw),FUN = sum,na.rm = T)

	df <- cbind(df1,df2[,2],.get.div(s))
	names(df) <- c("well","meanfiringrate","meanfiringrate_per_well","div")
	rownames(df)<- t(df["well"])
	df
}

.plot.isis.by.electrode<-function(s) {
	wells <- unique(s$cw)
	if (length(wells)>0) {
		for (well in wells) {
			active.electrodes <- which(s$cw == well & as.vector(unlist(lapply(s$isis,length))) >0)
			if (length(active.electrodes) > 0) {
				df <- list()
				for (i in 1:length(active.electrodes)) {
					df[[i]] <- cbind(s$isis[[active.electrodes[i]]],
					names(s$isis)[active.electrodes[i]])
				}
				df = do.call("rbind",df)
				colnames(df) <- c("isis","electrode")
				plateinfo <- plateinfo(s$layout$array)
				d1 <- expand.grid(col=1:plateinfo$n.elec.c,row=1:plateinfo$n.elec.r)
				all.electrodes <- sort(paste(well,"_", d1[,"row"],d1[,"col"],sep=""))
				layout.electrodes <- c(plateinfo$n.elec.r, plateinfo$n.elec.c)
				df <- data.frame(df)
				df$isis <- as.numeric(as.vector(df$isis))
				df$electrode <- as.character(as.vector(df$electrode))


				p1 <- histogram(~ isis | factor(electrode,levels=all.electrodes),
					 data = df, breaks = 10,
				 	main = paste("ISIs histogram plot for ", well, sep = "" ),
					layout = layout.electrodes,
					drop.unused.levels = FALSE)
				print(p1)
			
				p2 <- histogram(~ log(isis) | factor(electrode,levels=all.electrodes),
				 	data = df, breaks = 10,
				 	main = paste("log(ISIs) histogram plot for ", well, sep = "" ),
					layout = layout.electrodes,
					drop.unused.levels = FALSE)
				print(p2)
			}
		}
		p2
	}
	
}

.plot.mean.firingrate.by.electrode<-function(s) {
	wells <- unique(s$cw)
	if (length(wells)>0) {
		for (well in wells) {
			active.electrodes <- which(s$cw == well)
			df <- list()
			for (i in active.electrodes) {
				df[[i]] <- cbind(s$rates$times,
				s$rates$rates[,i],
				names(s$nspikes)[i])
			}
			df = do.call("rbind",df)
			maxy <- max(df[,2])
			colnames(df) <- c("time","meanfiringrate","electrode")
			plateinfo <- plateinfo(s$layout$array)
			d1 <- expand.grid(col=1:plateinfo$n.elec.c,row=1:plateinfo$n.elec.r)
			all.electrodes <- sort(paste(well,"_", d1[,"row"],d1[,"col"],sep=""))
			layout.electrodes <- c(plateinfo$n.elec.r, plateinfo$n.elec.c)
			df <- data.frame(df)
			
			p1 <- xyplot(meanfiringrate ~ time | factor(electrode,levels=all.electrodes),
				 data = df, 
				 main = paste("Mean Firing Rate per Second for Well ", well, ". Maximum firing rate:",maxy," Hz", sep = "" ),
				layout = layout.electrodes,type = "h",
				scales=list(
          				x=list(draw = FALSE),
          				y=list(draw = FALSE)),    
				drop.unused.levels = FALSE)
			print(p1)
		}
		p1
	}
	
}

#############################################
IGM.plot.mean.firingrate.by.eletrode.by.div<- function(s) {
	electrode.stats <- lapply(s, function(d){cbind(d$meanfiringrate,d$cw, .get.div(d))})
	electrode.stats.all <- do.call("rbind",electrode.stats)
	electrode.names <- row.names(electrode.stats.all)
	electrode.stats.all <- suppressWarnings(data.frame(cbind(electrode.names, electrode.stats.all[,1:3])))
	names(electrode.stats.all) <- c("electrode","meanfiringrate","well","div")
	electrode.stats.all$div <- as.numeric(as.vector(electrode.stats.all$div))
	electrode.stats.all$meanfiringrate <- as.numeric(as.vector(electrode.stats.all$meanfiringrate))
	electrode.stats.all$electrode <- as.character(as.vector(electrode.stats.all$electrode))


	wells <- unique(electrode.stats.all$well)
	if (length(wells)>0) {
		for (active.well in wells) {
			df <- electrode.stats.all[which(electrode.stats.all$well == active.well),]
			layout.info <- .get.electrode.layout(s[[1]],active.well)
			maxy <- max(df$meanfiringrate)

			p1 <- xyplot(meanfiringrate ~ div | factor(electrode,levels=layout.info$electrodes),
				 data = df, 
				 main = paste("Mean Firing Rate across DIV's for ", active.well, ". Maximum firing rate:",round(maxy,2)," Hz", sep = "" ),
				layout = layout.info$layout,
				#type = "h",
				#scales=list(
          			#	y=list(draw = FALSE)),    
				drop.unused.levels = FALSE)
			print(p1)
		}
	}
}

IGM.plot.mean.firingrate.by.well.by.div<- function(s) {
	well.stats <- lapply(s, function(d){d$well.stats})
	well.stats.all <- do.call("rbind",well.stats)
	plateinfo <- plateinfo(s[[1]]$layout$array)
	wells <- plateinfo$wells
	names(wells) <- wells #keep the names valid.
	wells.layout <- plateinfo$layout

	p1 <- xyplot(meanfiringrate ~ div | factor(well,levels=wells), data = well.stats.all, 
		main = "Mean Firing Rate across DIV's (Hz/electrode)",layout = wells.layout,
		drop.unused.levels = FALSE)
	print(p1)
	p2 <- xyplot(meanfiringrate_per_well ~ div | factor(well,levels=wells), data = well.stats.all, 
		main = "Mean Firing Rate across DIV's (Hz/well)",layout = wells.layout,
		drop.unused.levels = FALSE)
	print(p2)
	#return(list(p1=p1,p2=p2))
}

IGM.plot.plate.summary.for.spikes<-function(s,outputdir) {
	for (i in 1:length(s)) {
		basename <- get.file.basename(s[[i]]$file)
		spikePlotPath = paste(outputdir,"/",basename,"_spike_plot.pdf",sep="")
		pdf(file=spikePlotPath) 

		#layout 
  		p<-.plot.mealayout.1(s[[i]]$layout, use.names=T, cex=0.25)
  		title(main= paste( paste("Electrode Layout"), 
                     paste("file= ",  strsplit(basename(s[[i]]$file),".RData")[[1]][1], sep=""),sep='\n'))
		#MFR
		p<- .plot.meanfiringrate(s[[i]], main = "Mean Firing Rate by Plate (Hz)")
		#p<- plot(s[[i]],main = "Raster plots by channel", label.cells = FALSE, use.names = FALSE)
		p<- .plot.isis.by.plate(s[[i]])
  		p<-.channel.plot.by.well(s[[i]],resp="meanfiringrate", resp.label="Mean Firing Rate (Hz)")  
		p<-.plot.mean.firingrate.by.electrode(s[[i]])
		p<-.plot.isis.by.electrode(s[[i]])
		dev.off()
	}
}

write.plate.summary.for.spikes<-function(s,outputdir) {
	csvwell <- paste(outputdir,"/",get.project.plate.name(s[[1]]$file),"_well_spikes.csv",sep="")
	
	for (i in 1:length(s)) {
		div <- .get.div(s[[i]])
		basename <- get.file.basename(s[[i]]$file)
		csvfile <- paste(outputdir,"/",basename,"_spikes.csv",sep="")
		df <- .spike.summary.by.electrode(s[[i]])
		df2 <- .spike.summary.by.well(s[[i]])
		
		#recording time
		write.table(paste("recording time (s): [", paste(s[[i]]$rec.time[1],round(s[[i]]$rec.time[2]), sep=" ,"),
            	"]",sep=""),csvfile, sep=",", append=FALSE,row.names=FALSE,col.names=FALSE)
		write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
		write.table("Spike statistics for wells",	csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE) 
		df2 = cbind(rownames(df2),df2)
		suppressWarnings(write.table(df2,
                  csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=TRUE))
		suppressWarnings(write.table(cbind(df2,div),
                  csvwell, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE))

		write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)

		write.table("Spike statistics for electrodes",
      		csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE) 		
  		df = cbind(rownames(df),df)
		colnames(df)[1] <- "electrode"
      	#summary write data
      	suppressWarnings(write.table(df,
                  csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=TRUE))
	}
}

