on.off <-
function(directory=directory,id,visit,name.of.log.on.off,data)
{
		
	on.off.log <- read.csv(paste(directory,name.of.log.on.off,".csv",sep=""))
	on.off.log$id <- as.character(on.off.log$id)
	on.off.log <- on.off.log[on.off.log$id==id& on.off.log$visit==visit,]
	
	on.off.log$date.on <- paste(on.off.log$date.on.month,on.off.log$date.on.day,on.off.log$date.on.year,sep="/")
	on.off.log$time.on <- paste(on.off.log$time.on.hour,on.off.log$time.on.minute,on.off.log$time.on.seconds,sep=":")
	
	on.off.log$date.off <- paste(on.off.log$date.off.month,on.off.log$date.off.day,on.off.log$date.off.year,sep="/")
	on.off.log$time.off <- paste(on.off.log$time.off.hour,on.off.log$time.off.minute,on.off.log$time.off.seconds,sep=":")

	on.off.log$date.time.on <- paste(on.off.log$date.on, on.off.log$time.on, sep=" ")
	on.off.log$date.time.off <- paste(on.off.log$date.off, on.off.log$time.off, sep=" ")
	
	on.off.log$date.time.on <- strptime(on.off.log$date.time.on,"%m/%d/%Y %H:%M:%S")
	on.off.log$date.time.off <- strptime(on.off.log$date.time.off,"%m/%d/%Y %H:%M:%S")
	
	on.off.log$hours.on <- as.vector(difftime(strptime(on.off.log$date.time.off,format="%Y-%m-%d %H:%M:%S"),strptime(on.off.log$date.time.on,format="%Y-%m-%d %H:%M:%S"), units="hours"))
	
	#	if on/off times recorded - loop through and label time monitor is not worn
	if(dim(on.off.log)[1]>0)
	{
	data$off <- 1	
	for (t in (1:dim(on.off.log)[1]))
		{
	on <- strptime(on.off.log$date.time.on[t],"%Y-%m-%d %H:%M:%S")
    off <- strptime(on.off.log$date.time.off[t],"%Y-%m-%d %H:%M:%S")
		n <- dim(data)[1]
		inds <- (1:n)[((data$time>=on)&(data$time<=off))]
		if (length(inds)>0)
			data$off[inds] <- 0
			}
			if(dim(on.off.log)[1]==0)
			data$off <- "No.On.Off.Log"	
				}	#end bed loop
	return(data)
	}

