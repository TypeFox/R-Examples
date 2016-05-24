has.network.spikes <-function(nspikes) {
	return(!is.null(.as.data.frame.network.spikes(nspikes)))
}

IGM.summary.network.spikes<-function(e,nspikes,ns.E,sur){
	wells <- nspikes$wells
	names(wells) <- wells #keep the names valid.
	for (i in 1:length(wells)) {
		well <- wells[[i]]
		data <- nspikes$ns.all[[i]]$measures
		if (!is.null(data)) {
			indexes <- .names.to.indexes(names(e$spikes), well, allow.na = TRUE)
			electrodes <- names(e$spikes)[indexes]
			en.map <- matrix(0,length(electrodes),dim(data)[1])
			rownames(en.map) <- electrodes
			colnames(en.map) <-  as.character(data[,1])
			#beg <- floor(min(unlist(e$spikes[indexes])))
			for (ns in 1:dim(data)[1]) {
				current.ns <- data[ns,1]
				en.map[,ns] <- unlist(lapply(e$spikes[indexes],function(x) {length(x[x>current.ns & x<= current.ns + (nspikes$ns.all)[[i]]$ns.T])}))
			}

			en.map[en.map < ns.E] <- 0
			filtered.indexes <- which(colSums(en.map>=ns.E) >= (nspikes$ns.all)[[i]]$ns.N)

			en.map <- en.map[,filtered.indexes]
			if (length(filtered.indexes) == 1) {  #deal with R vector matrix problem
				dim(en.map) <- c(length(electrodes),length(filtered.indexes))
				colnames(en.map) <-  names(filtered.indexes)
				rownames(en.map) <- electrodes
			}

			

			if (dim(en.map)[2] > 0) {
				p <- data[filtered.indexes,1:2]
				dim(p) <- c(length(filtered.indexes),2)
				colnames(p)<-c("time","index")
				if (is.na(p[1,1])) {
					p<- NULL
				}
			} else {
				p <- NULL
			}
		} else {
			p <- NULL
		}
 
		ns <- list(counts = nspikes$ns.all[[i]]$counts, ns.N = (nspikes$ns.all)[[i]]$ns.N, ns.T = (nspikes$ns.all)[[i]]$ns.T)
	      class(ns) <- "ns"
        	m <- .mean.ns(ns, p, plot = FALSE, nrow = 4, ncol = 4, ask = FALSE, sur = sur)
        	if (is.null(m)) {
            	ns$brief <- c(n = 0, peak.m = NA, peak.sd = NA, durn.m = NA,durn.sd = NA,
			percent.of.spikes.in.ns = NA,
			mean.spikes.in.ns = NA,
			mean.insi = NA)
			ns$en.map <- NULL
			ns$brief.electrode <- NULL
        	} else {
            	ns$mean <- m$ns.mean
            	ns$measures <- m$measures
            	peak.val <- ns$measures[, "peak.val"]
            	durn <- ns$measures[, "durn"]
			insis <- diff(ns$measures[, "time"])
			if (length(insis) > 0) {
				mean.insis <- mean(insis)
			} else {
				mean.insis <- NA
			}
            	ns$brief <- c(n = nrow(ns$measures), peak.m = mean(peak.val), 
                	peak.sd = sd(peak.val), durn.m = mean(durn, na.rm = TRUE), 
                	durn.sd = sd(durn, na.rm = TRUE),
			percent.of.spikes.in.ns = 100*sum(en.map)/sum(e$nspikes[indexes]),
			mean.spikes.in.ns = sum(en.map) / nrow(ns$measures),
			mean.insis = mean.insis
			)
			if ( dim(en.map)[2] != dim(ns$measures)[1]) {
				en.map <- en.map[,as.character(ns$measures[,1])]
			}
			ns$en.map <- en.map

			#now briefs at electrode level
			features <- c("spikes","ns","spikes.in.ns","percent.of.spikes.in.ns",
				"mean.spikes.per.ns","sd.spikes.per.ns","mean.insis")
			en.brief <- matrix(0,dim(en.map)[1],length(features))
			rownames(en.brief) <- rownames(en.map)
			colnames(en.brief) <- features
			en.brief[,"spikes"] <- e$nspikes[indexes]
			en.brief[,"ns"] <-  rowSums(en.map>0)
			en.brief[,"spikes.in.ns"] <-  rowSums(en.map)
			en.brief[,"percent.of.spikes.in.ns"] <-  100 * en.brief[,"spikes.in.ns"] / en.brief[,"spikes"]

			en.brief[,"mean.spikes.per.ns"] <-  en.brief[,"spikes.in.ns"] / en.brief[,"ns"]
			en.brief[,"sd.spikes.per.ns"] <- unlist(lapply(rownames(en.brief), function(e) {
				temp <- sd(en.map[e,which(en.map[e,] >0)])
				temp[is.na(temp)] <- NaN
				temp
			}))


			en.brief[,"mean.insis"] <- unlist(lapply(rownames(en.brief), function(e) {
				insis <- diff(as.numeric(names(which(en.map[e,] >0))))
				m <- mean(insis)
				m
			}))
			
			en.brief[is.nan(en.brief)] <- NA

			ns$en.brief <- en.brief
        	}

		nspikes$ns.all[[i]] <- ns
	}
	nspikes
}

calculate.network.spikes <- function(e,sur=100, ns.N, ns.T) {
	#get well information
	plateinfo <- plateinfo(e$layout$array)
	wells <- plateinfo$wells
	names(wells) <- wells #keep the names valid.
	wells.layout <- plateinfo$layout
	##lets just use lapply
	ns.all<-lapply(wells, function(well) {
      	.compute.ns(e, ns.T = ns.T, ns.N = ns.N, sur=sur, whichcells = well)})
	return(list(wells = wells,ns.all = ns.all,wells.layout = wells.layout))
}

.as.data.frame.network.spikes <- function(nspikes) {
	#useful function
	mean.ns.to.xy <- function(ns, well) {
  		## Convert the mean network spike into a dataframe suitable for lattice
  		## graphics.
  		if (!is.null(ns$mean)) {
    			m = ns$mean
    			t = as.vector(time(m))
    			d <- data.frame(t = t, y = as.vector(m), well = well)
    			d
  		}
	}
	
	#get data in right format for lattice
	xys <- mapply(SIMPLIFY = FALSE, function(ns, well) {
  		mean.ns.to.xy(ns, well)
	}, nspikes$ns.all, names(nspikes$ns.all))

	d <- do.call("rbind", xys) #merge all the data frames.
	d	
}
IGM.xyplot.network.spikes <- function(nspikes) {
	p1 <- NULL
	## Produce the title for each well; note difference between '0' for zero
	## network spikes found in a well vs 'NA' (no spikes on that well).
	strip.names <- sapply(nspikes$wells, function(well) {
		n = nspikes$ns.all[[well]]
  		n.ns = n$brief[[1]]
  		paste(well, n.ns)
	})
	df <- .as.data.frame.network.spikes(nspikes)
	if (!is.null(df)) {
		p1 <- xyplot(y ~ 10*t | factor(well, levels = nspikes$wells), data = df, strip = strip.custom(factor.levels = strip.names),
      		main = "Mean NS",xlab = "Time (ms)", ylab = "# electrodes", type = "l", drop.unused.levels = FALSE,
       		layout = nspikes$wells.layout)
		print(p1)
	p1
	}
}
.active.wells.network.spikes<- function(nspikes) {
	active.wells <- nspikes
	active.wells$ns.all <- nspikes$ns.all[sort(nspikes$wells[sapply(nspikes$wells,function(well) {
		!is.na(nspikes$ns.all[[well]]$brief[[1]]) & nspikes$ns.all[[well]]$brief[[1]] > 0})
	])]
	active.wells
}

IGM.plot.active.wells.network.spikes<- function(nspikes) {
	active.wells <- .active.wells.network.spikes(nspikes)$ns.all
	if (length(active.wells)>0) {
		for (j in 1:length(active.wells)) {
			IGM.plot.ns(active.wells[[j]], main = names(active.wells)[j],ylab='Count', xlab='Time (s)')
			y<- as.vector(active.wells[[j]]$mean)
			plot(ts(y,start = c(-(length(y)-1)/2,1)), xlab='Time (ms)', ylab='Count', main=paste('Mean NS for', 
				names(active.wells)[j],sep = " "))
		}
	}
}

write.csv.network.spikes <-function(s,nspikes,outputdir) {
  csvwell <- paste(outputdir,"/",get.project.plate.name(s$file),"_ns.csv",sep="")
  div <- .get.div(s)
  active.wells <- .active.wells.network.spikes(nspikes)$ns.all
  if (length(active.wells) >0) {
    # sahar 10292014 - add genotype column and change newcol from 2 to 3
    newcol <-3
    #2 to peak.min and peak.max
    p <- length(active.wells[[1]]$brief) + length(active.wells[[1]]$mean) + newcol
    nsdata <- matrix(0,length(s$well),p)
    temp<-c() # Diana 10/2014
    # Diana Hall 10-31-2014 change
    length.temp.mean<-length(active.wells[[1]]$mean)
    for (j in 1:length(s$well)) {
      cur.well<-s$well[j]
      if ( is.element(cur.well, names(active.wells) ) ){
        temp <- active.wells[[cur.well]]
        nsdata[j,1:length(temp$brief)] <- temp$brief
        nsdata[j,length(temp$brief)+1] <- min(temp$measures[,"peak.val"])
        nsdata[j,length(temp$brief)+2] <- max(temp$measures[,"peak.val"])
        nsdata[j,length(temp$brief)+3] <- s$treatment[cur.well]
        nsdata[j,(length(temp$brief)+newcol+1):p] <- as.double(temp$mean)
        
      } else{
        temp$brief<-c(0, rep(NA,4),0, NA, NA)
        nsdata[j,1:length(temp$brief)] <- temp$brief
        nsdata[j,length(temp$brief)+1] <- NA
        nsdata[j,length(temp$brief)+2] <- NA
        #nsdata[j,length(temp$brief)+3] <- as.character( s$treatment[cur.well] )
        nsdata[j,length(temp$brief)+3] <- s$treatment[cur.well] 
        nsdata[j,(length(temp$brief)+newcol+1):p] <- rep(0, length.temp.mean)
      }
    }
    
    nsdata <- data.frame(nsdata)	
    names(nsdata)[1:length(temp$brief)] <- names(active.wells[[1]]$brief)
    names(nsdata)[(length(temp$brief)+1):(length(temp$brief)+newcol)] <- c("peak.min","peak.max","treatment")
    
    for (j in 1:(p-length(temp$brief)-newcol)) {
      names(nsdata)[j+newcol+length(temp$brief)] = paste("t",j,sep = '')
    }
    nsdata<- cbind(s$well,nsdata)
    names(nsdata)[1] <- "well"
    
    basename <- get.file.basename(s$file)
    csvfile= paste(outputdir,"/",basename,"_ns.csv",sep="")
    
    write.table(paste("file= ",  strsplit( basename(s$file),".RData")[[1]][1], sep=""),
                csvfile, sep=",", append=FALSE,row.names=FALSE,col.names=FALSE) 
    write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
    #recording time
    write.table(paste("recording time (s): [", paste(s$rec.time[1],round(s$rec.time[2]), sep=" ,"),
                      "]",sep=""),csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
    
    write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
    
    write.table("Network Spike analysis at well level",
                csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)   	
    suppressWarnings(write.table(nsdata,
                                 csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=TRUE))
    suppressWarnings(write.table(cbind(div,nsdata),
                                 csvwell, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE))
    
    
    write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
    write.table("Network Spike analysis at electrode level",
                csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)  
    
    
    en.df <- do.call("rbind",lapply(nspikes$ns.all, function(well) {
      temp <- well$en.brief
      temp
    }))
    en.df <- en.df[order(rownames(en.df)), ] 
    en.df <- cbind(rownames(en.df),en.df)
    colnames(en.df)[1] <- "electrode"
    suppressWarnings(write.table(en.df,
                                 csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=TRUE))
    
  }
}
