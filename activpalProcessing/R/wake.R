wake <-
function(directory=directory,id,visit,name.of.log.bed,data)
{
	bed.log <- read.csv(paste(directory,name.of.log.bed,".csv",sep=""))
	bed.log$id <- as.character(bed.log$id)
	bed.log <- bed.log[bed.log$id==id&bed.log$visit==visit,]
	
	
	bed.log$date.out <- paste(bed.log$date.out.month,bed.log$date.out.day,bed.log$date.out.year,sep="/")
	bed.log$time.out <- paste(bed.log$time.out.hour,bed.log$time.out.minute,bed.log$time.out.seconds,sep=":")

	bed.log$date.in <- paste(bed.log$date.in.month,bed.log$date.in.day,bed.log$date.in.year,sep="/")
	bed.log$time.in <- paste(bed.log$time.in.hour,bed.log$time.in.minute,bed.log$time.in.seconds,sep=":")

	bed.log$date.time.out <- paste(bed.log$date.out, bed.log$time.out, sep=" ")
	bed.log$date.time.in <- paste(bed.log$date.in, bed.log$time.in, sep=" ")
	
	bed.log$date.time.out <- strptime(bed.log$date.time.out,"%m/%d/%Y %H:%M:%S")
	bed.log$date.time.in <- strptime(bed.log$date.time.in,"%m/%d/%Y %H:%M:%S")
	
	bed.log$hours.up <- as.vector(difftime(strptime(bed.log$date.time.in,format="%Y-%m-%d %H:%M:%S"),strptime(bed.log$date.time.out,format="%Y-%m-%d %H:%M:%S"), units="hours"))

#	if bed times recorded - loop through and label time in bed
	if(dim(bed.log)[1]>0)
	{
	data$in.bed <- 1
	data$day.for.wearer <- NA	
	bed.log$day.for.wearer <- 1:dim(bed.log)[1]
		#	t <- 1
	for (t in (1:dim(bed.log)[1]))
		{
	wake <- strptime(bed.log$date.time.out[t],"%Y-%m-%d %H:%M:%S")
    bed <- strptime(bed.log$date.time.in[t],"%Y-%m-%d %H:%M:%S")
		n <- dim(data)[1]
		inds <- (1:n)[((data$time>=wake)&(data$time<=bed))]

		if (length(inds)>0)
			data$in.bed[inds] <- 0
			data$day.for.wearer[inds] <- bed.log$day.for.wearer[t]
			}
			if(dim(bed.log)[1]==0)
			data$in.bed <- "No.Bed.Log"	
				}	#end bed loop
	
	return(data)
	}

