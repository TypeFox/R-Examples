process.AP <-
function(directory,name.of.log.subjects,name.of.log.bed=NULL,name.of.log.on.off=NULL)
{
	counter <- 1
	
	subs <- identifySubjects(directory,name.of.log.subjects)
	sn <- length(subs)
	visit <- identifyVisits(directory,name.of.log.subjects)
	vn <- length(visit)
	study <- identifyStudy(directory,name.of.log.subjects)
					
 for (s in subs)
	{
	for (v in visit)
		{
			print(paste("Study ID",s,sep=" "))
			print(v)
			print(study)
			print(paste(counter,"of",sn*vn,sep=" "))
					
	 temp <- Sys.glob(paste(directory,study,"_",s,"_",v,".csv",sep=""))
		
		if (length(temp)>0)
			{
		file.name.and.path <- paste(directory,study,"_",s,"_",v,".csv",sep="")
		
		temp <- activpal.file.reader(file.name.and.path)
		
#		head(temp)
#		dim(temp)

		n <- dim(temp)[1]		
	
		temp <- second.by.second(temp)

#	loop to label sleep time
	
	bed.log.temp <- Sys.glob(paste(directory,name.of.log.bed,".csv",sep=""))
	
	temp$in.bed <- 1
	temp$day.for.wearer <- NA	
	
	if(length(bed.log.temp>0))
{ 
	bed.log <- read.csv(bed.log.temp)
	bed.log$id <- as.character(bed.log$id)
	bed.log <- bed.log[bed.log$id==s&bed.log$visit==v,]
	
	if(dim(bed.log)[1]>0)
	{
		
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

	bed.log$day.for.wearer <- 1:dim(bed.log)[1]

	for (t in (1:dim(bed.log)[1]))
		{
	wake <- strptime(bed.log$date.time.out[t],"%Y-%m-%d %H:%M:%S")
    bed <- strptime(bed.log$date.time.in[t],"%Y-%m-%d %H:%M:%S")
		n <- dim(temp)[1]
		inds <- (1:n)[((temp$time>=wake)&(temp$time<=bed))]

		if (length(inds)>0)
			temp$in.bed[inds] <- 0
			temp$day.for.wearer[inds] <- bed.log$day.for.wearer[t]
			}

			inds.time.awake <- (1:(dim(temp)[1]))[temp$in.bed==0]
			a <- length(inds.time.awake)
			if(a==0)
			temp$in.bed <- "AP and bed.log do not match"
				
					}
		
				if(dim(bed.log)[1]==0)
				temp$in.bed <- "Subject/Visit not in bed.log"					
					}	#end sleep time loop
					
			if(length(bed.log.temp)==0)
			temp$in.bed <- "No.Bed.Log"	
				
				
#	loop to remove label on/off time
	on.off.log.temp <- Sys.glob(paste(directory,name.of.log.on.off,".csv",sep=""))
	
	temp$off <- 1	

	if(length(on.off.log.temp>0))
{ 
	on.off.log <- read.csv(on.off.log.temp)
	on.off.log$id <- as.character(on.off.log$id)
	on.off.log <- on.off.log[on.off.log$id==s& on.off.log$visit==v,]
	
	if(dim(on.off.log)[1]>0)
	{
		
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
	
	for (t in (1:dim(on.off.log)[1]))
		{
	on <- strptime(on.off.log$date.time.on[t],"%Y-%m-%d %H:%M:%S")
    off <- strptime(on.off.log$date.time.off[t],"%Y-%m-%d %H:%M:%S")
		n <- dim(temp)[1]
		inds <- (1:n)[((temp$time>=on)&(temp$time<=off))]
		if (length(inds)>0)
			temp$off[inds] <- 0
			}
			
			inds.time.worn <- (1:(dim(temp)[1]))[temp$off==0]
			i <- length(inds.time.worn)
			if(i==0)
			temp$off <- "AP and on.off.log do not match"
			
			}	
				
			if(dim(on.off.log)[1]==0)
		temp$off <- "Subject/Visit not in on.off.log"		
					}	#end on/off loop

			if(length(on.off.log.temp)==0)
			temp$off <- "No.On.Off.Log"		
							
			temp$counter <- 1
			
#	get step count - cumulative steps reported in file so 	need to figure out total steps/day
			
	d <- dim(temp)[1]
		
	steps.inds <- c((1:d)[temp$date[-1]!=temp$date[-d]],d)
	steps.1 <- temp[steps.inds,]
	steps.2 <- as.vector(steps.1$steps)
	d.2 <- length(steps.2)
	steps.3 <- data.frame(date=temp$date[steps.inds],steps=c(steps.2[1],steps.2[-1]-steps.2[-d.2]))
	
	####		
		
	if(is.numeric(temp$in.bed)==T&is.numeric(temp$off)==T)
	{
		inds.time.worn <- (1:(dim(temp)[1]))[temp$off==0]
		inds.first.time.worn <- inds.time.worn[1]
		first.time.worn <- temp$time[inds.first.time.worn]
		
		inds.last.time.worn <- inds.time.worn[i]
		last.time.worn <- temp$time[inds.last.time.worn]
		
		only.measurement.period <- temp[temp$time>=first.time.worn&temp$time<=last.time.worn,]
		dim(only.measurement.period)
		only.measurement.period.awake <- only.measurement.period[only.measurement.period$in.bed==0,]
		dim(only.measurement.period.awake)
		
		sleep.wake.wear.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date=unique(only.measurement.period$date),
		awake.hours=tapply(only.measurement.period$in.bed==0, only.measurement.period$date,sum)/3600,
		total.sleep.hours=tapply(only.measurement.period$in.bed==1, only.measurement.period$date,sum)/3600,
		total.wear.hours=tapply(only.measurement.period$off==0, only.measurement.period$date,sum)/3600,
		non.wear.hours=tapply(only.measurement.period$off==1, only.measurement.period$date,sum)/3600,
		hours.awake.worn=(tapply(only.measurement.period$in.bed==0&only.measurement.period$off==0,only.measurement.period$date,sum))/3600,
		hours.awake.not.worn=(tapply(only.measurement.period$in.bed==0&only.measurement.period$off==1,only.measurement.period$date,sum))/3600,
		hours.sleep.worn=(tapply(only.measurement.period$in.bed==1&only.measurement.period$off==0,only.measurement.period$date,sum))/3600,
		hours.sleep.not.worn=(tapply(only.measurement.period$in.bed==1&only.measurement.period$off==1,only.measurement.period$date,sum))/3600
		)		
		
		# make file with bed and off time cleaned out
		data <- temp[temp$in.bed==0&temp$off==0,]
		
		if(dim(data)[1]>1)
	{
#		head(data)	

#	get step count - cumulative steps reported in file so 		need to figure out total steps/day

	date <- unique(data$date)
	dd <- dim(steps.3)[1]
	
	inds <- vector(length=0)
	
	for(i in (1:length(date)))
	{
	inds <- c(inds,(1:dd)[date[i]==steps.3$date])
		}
	
	steps <- steps.3$steps[inds]
	
		# make .csv file with PA and SB variables per day
		
		results.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date=unique(data$date),
	
	hours.awake.worn = tapply(data$off==0,data$date,sum)/3600,
		
	met.hours = tapply(data$met.hours,data$date,sum),	step.count = steps,

	sed.mins = tapply(data$ap.posture,data$date,sed.min.AP),
	stand.mins = tapply(data$ap.posture,data$date,stand.min.AP),
	step.mins = tapply(data$ap.posture,data$date,step.min.AP),

	lit.mins = tapply(data$mets,data$date,lit.min.AP),
	mvpa.mins = tapply(data$mets,data$date,mvpa.min.AP),
	
	breaks = tapply(data$ap.posture,data$date,breaks.AP),
	break.rate = tapply(data$ap.posture,data$date,breaks.AP)/(tapply(data$ap.posture,data$date,sed.min.AP)),
	
	guideline.minutes = tapply(data$one.min.mets,data$date,guideline.bouts.min),
	num.guideline.bouts = tapply(data$one.min.mets,data$date,guideline.bouts.num),

	min.in.sed.30 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.min,n=30),
	min.in.sed.60 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.min,n=60),
	
	num.bouts.in.sed.30 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.num,n=30),
	num.bouts.in.sed.60 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.num,n=60))
	
	}
	
	if(dim(data)[1]==0)
	{

		results.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date="No Valid Wear Hours.  Check that AP file matches with on.off.log and bed.log",
	
	hours.awake.worn = NA,
		
	met.hours = NA,
	step.count = NA,

	sed.mins = NA,
	stand.mins = NA,
	step.mins = NA,

	lit.mins = NA,
	mvpa.mins = NA,
	
	breaks = NA,
	break.rate = NA,
	
	guideline.minutes = NA,
	num.guideline.bouts = NA,

	min.in.sed.30 = NA,
	min.in.sed.60 = NA,
	
	num.bouts.in.sed.30 = NA,
	num.bouts.in.sed.60 = NA)
	
	}
		}

####

	if(is.numeric(temp$in.bed)==F&is.numeric(temp$off)==T)
	{
		inds.time.worn <- (1:(dim(temp)[1]))[temp$off==0]
		i <- length(inds.time.worn)
		inds.first.time.worn <- inds.time.worn[1]
		first.time.worn <- temp$time[inds.first.time.worn]
		
		inds.last.time.worn <- inds.time.worn[i]
		last.time.worn <- temp$time[inds.last.time.worn]
		
		only.measurement.period <- temp[temp$time>=first.time.worn&temp$time<=last.time.worn,]
		dim(only.measurement.period)
		
		sleep.wake.wear.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date=unique(only.measurement.period$date),
		awake.hours="No valid bed.log data",
		sleep.hours="No valid bed.log data",
		total.wear.hours=tapply(only.measurement.period$off==0, only.measurement.period$date,sum)/3600,
		total.non.wear.hours=tapply(only.measurement.period$off==1, only.measurement.period$date,sum)/3600,
		hours.awake.worn="No valid bed.log data",
		hours.awake.not.worn="No valid bed.log data",
		hours.sleep.worn="No valid bed.log data",
		hours.sleep.not.worn="No valid bed.log data"
		)		
		
		# make file with bed and off time cleaned out
		data <- temp[temp$off==0,]
		
		if(dim(data)[1]>1)
	{
#		head(data)	

	#	get step count - cumulative steps reported in file so 		need to figure out total steps/day

	date <- unique(data$date)
	dd <- dim(steps.3)[1]
	
	inds <- vector(length=0)
	
	for(i in (1:length(date)))
	{
	inds <- c(inds,(1:dd)[date[i]==steps.3$date])
		}
	
	steps <- steps.3$steps[inds]
	
		# make .csv file with PA and SB variables per day
		
		results.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date=unique(data$date),
	
	hours.awake.worn = tapply(data$off==0,data$date,sum)/3600,
	
	met.hours = tapply(data$met.hours,data$date,sum),	step.count = steps,
	
	sed.mins = tapply(data$ap.posture,data$date,sed.min.AP),
	stand.mins = tapply(data$ap.posture,data$date,stand.min.AP),
	step.mins = tapply(data$ap.posture,data$date,step.min.AP),

	lit.mins = tapply(data$mets,data$date,lit.min.AP),
	mvpa.mins = tapply(data$mets,data$date,mvpa.min.AP),
	
	breaks = tapply(data$ap.posture,data$date,breaks.AP),
	break.rate = tapply(data$ap.posture,data$date,breaks.AP)/(tapply(data$ap.posture,data$date,sed.min.AP)),
	
	guideline.minutes = tapply(data$one.min.mets,data$date,guideline.bouts.min),
	num.guideline.bouts = tapply(data$one.min.mets,data$date,guideline.bouts.num),

	min.in.sed.30 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.min,n=30),
	min.in.sed.60 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.min,n=60),
	
	num.bouts.in.sed.30 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.num,n=30),
	num.bouts.in.sed.60 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.num,n=60))
	
	}
	
	if(dim(data)[1]==0)
	{

	results.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date="No Valid Wear Hours.  Check that AP file matches with on.off.log and bed.log",
	
	hours.awake.worn = NA,
	
	met.hours = NA,
	step.count = NA,
	
	sed.mins = NA,
	stand.mins = NA,
	step.mins = NA,

	lit.mins = NA,
	mvpa.mins = NA,
	
	breaks = NA,
	break.rate = NA,
	
	guideline.minutes = NA,
	num.guideline.bouts = NA,

	min.in.sed.30 = NA,
	min.in.sed.60 = NA,
	
	num.bouts.in.sed.30 = NA,
	num.bouts.in.sed.60 = NA)
	
	}
		}
		
####

	if(is.numeric(temp$in.bed)==T&is.numeric(temp$off)==F)
	{
		only.measurement.period.awake <- temp[temp$in.bed==0,]
		dim(only.measurement.period.awake)
		
		sleep.wake.wear.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date=unique(only.measurement.period.awake$date),
		awake.hours=tapply(only.measurement.period.awake$in.bed==0, only.measurement.period.awake$date,sum)/3600,
		sleep.hours=tapply(only.measurement.period.awake$in.bed==1, only.measurement.period.awake$date,sum)/3600,
		total.wear.hours="No valid on.off.log data",
		total.non.wear.hours="No valid on.off.log data",
		hours.awake.worn="No valid on.off.log data",
		hours.awake.not.worn="No valid on.off.log data",
		hours.sleep.worn="No valid on.off.log data",
		hours.sleep.not.worn="No valid on.off.log data"
		)		
		
		# make file with bed and off time cleaned out
		data <- temp[temp$in.bed==0,]
		
		if(dim(data)[1]>1)
	{
#		head(data)	

		#	get step count - cumulative steps reported in file so 		need to figure out total steps/day

	date <- unique(data$date)
	dd <- dim(steps.3)[1]
	
	inds <- vector(length=0)
	
	for(i in (1:length(date)))
	{
	inds <- c(inds,(1:dd)[date[i]==steps.3$date])
		}
	
	steps <- steps.3$steps[inds]
	
	# make .csv file with PA and SB variables per day
		
	results.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date=unique(data$date),
	
	hours.awake.worn = tapply(data$off==0,data$date,sum)/3600,
		
	met.hours = tapply(data$met.hours,data$date,sum),	step.count = steps,

	sed.mins = tapply(data$ap.posture,data$date,sed.min.AP),
	stand.mins = tapply(data$ap.posture,data$date,stand.min.AP),
	step.mins = tapply(data$ap.posture,data$date,step.min.AP),

	lit.mins = tapply(data$mets,data$date,lit.min.AP),
	mvpa.mins = tapply(data$mets,data$date,mvpa.min.AP),
	
	breaks = tapply(data$ap.posture,data$date,breaks.AP),
	break.rate = tapply(data$ap.posture,data$date,breaks.AP)/(tapply(data$ap.posture,data$date,sed.min.AP)),
	
	guideline.minutes = tapply(data$one.min.mets,data$date,guideline.bouts.min),
	num.guideline.bouts = tapply(data$one.min.mets,data$date,guideline.bouts.num),

	min.in.sed.30 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.min,n=30),
	min.in.sed.60 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.min,n=60),
	
	num.bouts.in.sed.30 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.num,n=30),
	num.bouts.in.sed.60 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.num,n=60))
	
	}
	
	if(dim(data)[1]==0)
	{

		results.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date="No Valid Wear Hours.  Check that AP file matches with on.off.log and bed.log",
	
	hours.awake.worn = NA,
	
	met.hours = NA,
	step.count = NA,
	
	sed.mins = NA,
	stand.mins = NA,
	step.mins = NA,

	lit.mins = NA,
	mvpa.mins = NA,
	
	breaks = NA,
	break.rate = NA,
	
	guideline.minutes = NA,
	num.guideline.bouts = NA,

	min.in.sed.30 = NA,
	min.in.sed.60 = NA,
	
	num.bouts.in.sed.30 = NA,
	num.bouts.in.sed.60 = NA)
	
	}
		}		
		
####

	if(is.numeric(temp$in.bed)==F&is.numeric(temp$off)==F)
	{

		sleep.wake.wear.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date=unique(temp$date),
		awake.hours="No valid bed.log or on.off.log data",
		sleep.hours="No valid bed.log or on.off.log data",
		total.wear.hours="No valid bed.log or on.off.log data",
		total.non.wear.hours="No valid bed.log or on.off.log data",
		hours.awake.worn="No valid bed.log or on.off.log data",
		hours.awake.not.worn="No valid bed.log or on.off.log data",
		hours.sleep.worn="No valid bed.log or on.off.log data",
		hours.sleep.not.worn="No valid bed.log or on.off.log data"
		)				
		
		# make file with bed and off time cleaned out
		data <- temp
		if(dim(data)[1]>1)
	{
#		head(data)	

	#	get step count - cumulative steps reported in file so 		need to figure out total steps/day

	date <- unique(data$date)
	dd <- dim(steps.3)[1]
	
	inds <- vector(length=0)
	
	for(i in (1:length(date)))
	{
	inds <- c(inds,(1:dd)[date[i]==steps.3$date])
		}
	
	steps <- steps.3$steps[inds]
	
		# make .csv file with PA and SB variables per day
		
		results.table <- data.frame(study=study,sub=as.numeric(s),visit=v,date=unique(data$date),
	
	hours.awake.worn = tapply(data$off==0,data$date,sum)/3600,
		
	met.hours = tapply(data$met.hours,data$date,sum),	step.count = steps,

	sed.mins = tapply(data$ap.posture,data$date,sed.min.AP),
	stand.mins = tapply(data$ap.posture,data$date,stand.min.AP),
	step.mins = tapply(data$ap.posture,data$date,step.min.AP),

	lit.mins = tapply(data$mets,data$date,lit.min.AP),
	mvpa.mins = tapply(data$mets,data$date,mvpa.min.AP),
			
	breaks = tapply(data$ap.posture,data$date,breaks.AP),
	break.rate = tapply(data$ap.posture,data$date,breaks.AP)/(tapply(data$ap.posture,data$date,sed.min.AP)),
	
	guideline.minutes = tapply(data$one.min.mets,data$date,guideline.bouts.min),
	num.guideline.bouts = tapply(data$one.min.mets,data$date,guideline.bouts.num),

	min.in.sed.30 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.min,n=30),
	min.in.sed.60 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.min,n=60),
	
	num.bouts.in.sed.30 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.num,n=30),
	num.bouts.in.sed.60 = tapply(data$ap.posture,data$date,prolonged.sed.bouts.num,n=60))
	
	}
		}
	
	results.table$percent.of.hours.awake.worn.sed <- results.table$sed.mins/(results.table$hours.awake.worn*60)
	results.table$percent.of.hours.awake.worn.lit <- results.table$lit.mins/(results.table$hours.awake.worn*60)
	results.table$percent.of.hours.awake.worn.mvpa <- results.table$mvpa.mins/(results.table$hours.awake.worn*60)	
	
	rt <- dim(results.table)[1]
	inds.sed <- results.table$percent.of.hours.awake.worn.sed
	inds.inf.sed <- (1:rt)[inds.sed=="Inf"]
	results.table$percent.of.hours.awake.worn.sed[inds.inf.sed] <- NA

	inds.lit <- results.table$percent.of.hours.awake.worn.lit
	inds.inf.lit <- (1:rt)[inds.lit=="Inf"]
	results.table$percent.of.hours.awake.worn.lit[inds.inf.lit] <- NA

	inds.mvpa <- results.table$percent.of.hours.awake.worn.mvpa
	inds.inf.mvpa <- (1:rt)[inds.mvpa=="Inf"]
	results.table$percent.of.hours.awake.worn.mvpa[inds.inf.mvpa] <- NA
	
	means.table <- data.frame(study=study,sub=as.numeric(s),visit=v,
				
		hours.awake.worn = mean(results.table$hours.awake.worn,na.rm=T),
		sd.hours.awake.worn = sd(results.table$hours.awake.worn,na.rm=T),
		low.hours.awake.worn = mean(results.table$hours.awake.worn,na.rm=T)-1.96*(sd(results.table$hours.awake.worn,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.hours.awake.worn = mean(results.table$hours.awake.worn,na.rm=T)+1.96*(sd(results.table$hours.awake.worn,na.rm=T)/(sqrt(dim(results.table)[1]))),

		met.hours = mean(results.table$met.hours,na.rm=T),	
		sd.met.hours = sd(results.table$met.hours,na.rm=T),	
		low.met.hours = mean(results.table$met.hours,na.rm=T)-1.96*(sd(results.table$met.hours,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.met.hours = mean(results.table$met.hours,na.rm=T)+1.96*(sd(results.table$met.hours,na.rm=T)/(sqrt(dim(results.table)[1]))),

		step.count = mean(results.table$step.count,na.rm=T),
		sd.step.count = sd(results.table$step.count,na.rm=T),
		low.step.count = mean(results.table$step.count,na.rm=T)-1.96*(sd(results.table$step.count,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.step.count = mean(results.table$step.count,na.rm=T)+1.96*(sd(results.table$step.count,na.rm=T)/(sqrt(dim(results.table)[1]))),

		sed.mins = mean(results.table$sed.mins,na.rm=T),
		sd.sed.mins = sd(results.table$sed.mins,na.rm=T),
		low.sed.mins = mean(results.table$sed.mins,na.rm=T)-1.96*(sd(results.table$sed.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.sed.mins = mean(results.table$sed.mins,na.rm=T)+1.96*(sd(results.table$sed.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),

		stand.mins = mean(results.table$stand.mins,na.rm=T),
		sd.stand.mins = sd(results.table$stand.mins,na.rm=T),
		low.stand.mins = mean(results.table$stand.mins,na.rm=T)-1.96*(sd(results.table$stand.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.stand.mins = mean(results.table$stand.mins,na.rm=T)+1.96*(sd(results.table$stand.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),

		step.mins = mean(results.table$step.mins,na.rm=T),
		sd.step.mins = sd(results.table$step.mins,na.rm=T),
		low.step.mins = mean(results.table$step.mins,na.rm=T)-1.96*(sd(results.table$step.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.step.mins = mean(results.table$step.mins,na.rm=T)+1.96*(sd(results.table$step.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),

		lit.mins = mean(results.table$lit.mins,na.rm=T),
		sd.lit.mins = sd(results.table$lit.mins,na.rm=T),
		low.lit.mins = mean(results.table$lit.mins,na.rm=T)-1.96*(sd(results.table$lit.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.lit.mins = mean(results.table$lit.mins,na.rm=T)+1.96*(sd(results.table$lit.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),

		mvpa.mins = mean(results.table$mvpa.mins,na.rm=T),
		sd.mvpa.mins = sd(results.table$mvpa.mins,na.rm=T),
		low.mvpa.mins = mean(results.table$mvpa.mins,na.rm=T)-1.96*(sd(results.table$mvpa.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.mvpa.mins = mean(results.table$mvpa.mins,na.rm=T)+1.96*(sd(results.table$mvpa.mins,na.rm=T)/(sqrt(dim(results.table)[1]))),
			
		breaks = mean(results.table$breaks,na.rm=T),
		sd.breaks = sd(results.table$breaks,na.rm=T),
		low.breaks = mean(results.table$breaks,na.rm=T)-1.96*(sd(results.table$breaks,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.breaks = mean(results.table$breaks,na.rm=T)+1.96*(sd(results.table$breaks,na.rm=T)/(sqrt(dim(results.table)[1]))),

		break.rate = mean(results.table$break.rate,na.rm=T),
		sd.break.rate = sd(results.table$break.rate,na.rm=T),
		low.break.rate = mean(results.table$break.rate,na.rm=T)-1.96*(sd(results.table$break.rate,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.break.rate = mean(results.table$break.rate,na.rm=T)+1.96*(sd(results.table$break.rate,na.rm=T)/(sqrt(dim(results.table)[1]))),
	
		guideline.minutes = mean(results.table$guideline.minutes,na.rm=T),
		sd.guideline.minutes = sd(results.table$guideline.minutes,na.rm=T),
		low.guideline.minutes = mean(results.table$guideline.minutes,na.rm=T)-1.96*(sd(results.table$guideline.minutes,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.guideline.minutes = mean(results.table$guideline.minutes,na.rm=T)+1.96*(sd(results.table$guideline.minutes,na.rm=T)/(sqrt(dim(results.table)[1]))),

		num.guideline.bouts = mean(results.table$num.guideline.bouts,na.rm=T),
		sd.num.guideline.bouts = sd(results.table$num.guideline.bouts,na.rm=T),
		low.num.guideline.bouts = mean(results.table$num.guideline.bouts,na.rm=T)-1.96*(sd(results.table$num.guideline.bouts,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.num.guideline.bouts = mean(results.table$num.guideline.bouts,na.rm=T)+1.96*(sd(results.table$num.guideline.bouts,na.rm=T)/(sqrt(dim(results.table)[1]))),

		min.in.sed.30 = mean(results.table$min.in.sed.30,na.rm=T),
		sd.min.in.sed.30 = sd(results.table$min.in.sed.30,na.rm=T),
		low.min.in.sed.30 = mean(results.table$min.in.sed.30,na.rm=T)-1.96*(sd(results.table$min.in.sed.30,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.min.in.sed.30 = mean(results.table$min.in.sed.30,na.rm=T)+1.96*(sd(results.table$min.in.sed.30,na.rm=T)/(sqrt(dim(results.table)[1]))),

		min.in.sed.60 = mean(results.table$min.in.sed.60,na.rm=T),
		sd.min.in.sed.60 = sd(results.table$min.in.sed.60,na.rm=T),
		low.min.in.sed.60 = mean(results.table$min.in.sed.60,na.rm=T)-1.96*(sd(results.table$min.in.sed.60,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.min.in.sed.60 = mean(results.table$min.in.sed.60,na.rm=T)+1.96*(sd(results.table$min.in.sed.60,na.rm=T)/(sqrt(dim(results.table)[1]))),
	
		num.bouts.in.sed.30 = mean(results.table$num.bouts.in.sed.30,na.rm=T),
		sd.num.bouts.in.sed.30 = sd(results.table$num.bouts.in.sed.30,na.rm=T),
		low.num.bouts.in.sed.30 = mean(results.table$num.bouts.in.sed.30,na.rm=T)-1.96*(sd(results.table$num.bouts.in.sed.30,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.num.bouts.in.sed.30 = mean(results.table$num.bouts.in.sed.30,na.rm=T)+1.96*(sd(results.table$num.bouts.in.sed.30,na.rm=T)/(sqrt(dim(results.table)[1]))),

		num.bouts.in.sed.60 = mean(results.table$num.bouts.in.sed.60,na.rm=T),
		sd.num.bouts.in.sed.60 = sd(results.table$num.bouts.in.sed.60,na.rm=T),
		low.num.bouts.in.sed.60 = mean(results.table$num.bouts.in.sed.60,na.rm=T)-1.96*(sd(results.table$num.bouts.in.sed.60,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.num.bouts.in.sed.60 = mean(results.table$num.bouts.in.sed.60,na.rm=T)+1.96*(sd(results.table$num.bouts.in.sed.60,na.rm=T)/(sqrt(dim(results.table)[1]))),

		percent.of.hours.awake.worn.sed = mean(results.table$percent.of.hours.awake.worn.sed,na.rm=T),
		sd.percent.of.hours.awake.worn.sed = sd(results.table$percent.of.hours.awake.worn.sed,na.rm=T),
		low.percent.of.hours.awake.worn.sed = mean(results.table$percent.of.hours.awake.worn.sed,na.rm=T)-1.96*(sd(results.table$percent.of.hours.awake.worn.sed,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.percent.of.hours.awake.worn.sed = mean(results.table$percent.of.hours.awake.worn.sed,na.rm=T)+1.96*(sd(results.table$percent.of.hours.awake.worn.sed,na.rm=T)/(sqrt(dim(results.table)[1]))),

		percent.of.hours.awake.worn.lit = mean(results.table$percent.of.hours.awake.worn.lit,na.rm=T),
		sd.percent.of.hours.awake.worn.lit = sd(results.table$percent.of.hours.awake.worn.lit,na.rm=T),
		low.percent.of.hours.awake.worn.lit = mean(results.table$percent.of.hours.awake.worn.lit,na.rm=T)-1.96*(sd(results.table$percent.of.hours.awake.worn.lit,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.percent.of.hours.awake.worn.lit = mean(results.table$percent.of.hours.awake.worn.lit,na.rm=T)+1.96*(sd(results.table$percent.of.hours.awake.worn.lit,na.rm=T)/(sqrt(dim(results.table)[1]))),

		percent.of.hours.awake.worn.mvpa = mean(results.table$percent.of.hours.awake.worn.mvpa,na.rm=T),
		sd.percent.of.hours.awake.worn.mvpa = sd(results.table$percent.of.hours.awake.worn.mvpa,na.rm=T),
		low.percent.of.hours.awake.worn.mvpa = mean(results.table$percent.of.hours.awake.worn.mvpa,na.rm=T)-1.96*(sd(results.table$percent.of.hours.awake.worn.mvpa,na.rm=T)/(sqrt(dim(results.table)[1]))),
		up.percent.of.hours.awake.worn.mvpa = mean(results.table$percent.of.hours.awake.worn.mvpa,na.rm=T)+1.96*(sd(results.table$percent.of.hours.awake.worn.mvpa,na.rm=T)/(sqrt(dim(results.table)[1])))

			)

	if (counter==1)
	write.table(results.table,file=paste(directory,"results.table.csv",sep=""),sep=",",row.names=F,col.names=T,append=F)
	
	if (counter==1)
	
	write.table(sleep.wake.wear.table,file=paste(directory,"sleep.wake.wear.table.csv",sep=""),sep=",",row.names=F,col.names=T,append=F)

	if (counter==1)
	write.table(means.table,file=paste(directory,"means.table.csv",sep=""),sep=",",row.names=F,col.names=T,append=F)
	
	
	if (counter>1)
	  write.table(results.table,file=paste(directory,"results.table.csv",sep=""),sep=",",row.names=F,col.names=F,append=T)

	if (counter>1)
	write.table(sleep.wake.wear.table,file=paste(directory,"sleep.wake.wear.table.csv",sep=""),sep=",",row.names=F,col.names=F,append=T)
	  
	if (counter>1)
	write.table(means.table,file=paste(directory,"means.table.csv",sep=""),sep=",",row.names=F,col.names=F,append=T)
	
			
			
	counter <- counter+1
	
		}
			}
				}
				
					}

