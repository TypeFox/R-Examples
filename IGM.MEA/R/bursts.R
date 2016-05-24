IGM.plot.plate.summary.for.bursts<-function(s,outputdir,parameters) {
	for (i in (1:length(s))){
		basename <- get.file.basename(s[[i]]$file)
		burstPlotPath = paste(outputdir,"/",basename,"_burst_plot.pdf",sep="")
		pdf(file=burstPlotPath) 
		#layout 
  		p<-.plot.mealayout.1(s[[i]]$layout, use.names=T, cex=0.25)
  		title(main= paste( paste("Electrode Layout"), 
                     paste("file= ",  strsplit( basename(s[[i]]$file),".RData")[[1]][1], sep=""),                             
                     sep='\n'))

		# Add distribution plots 
		
		# Perform IBI distribution analysis
		if (parameters$burst.distribution.IBI$perform){
		  
		  feature="IBI"; print ("Running IBI distribution analysis.")
		  params=parameters$burst.distribution.IBI
		  p<-IGM.plot.distributions(s[[i]],minVals=params$min.cases,xlimit=params$x.lim,binsInSec=params$bins.in.seg,feature=feature,
		                        filterValuesByMin=params$filter.by.min,minValues=params$min.values,perWell=params$per.well,outputdir=outputdir)}
		
		# Perform ISI distribution analysis
		if (parameters$burst.distribution.ISI$perform){
		  
		  feature="ISI"; print ("Running ISI distribution analysis.")
		  params=parameters$burst.distribution.ISI
		  p<-IGM.plot.distributions(s[[i]],minVals=params$min.cases,xlimit=params$x.lim,binsInSec=params$bins.in.seg,feature=feature,
		                        filterValuesByMin=params$filter.by.min,minValues=params$min.values,perWell=params$per.well,outputdir=outputdir)}
		
		# Perform nSpikes distribution analysis
		if (parameters$burst.distribution.nSpikes$perform){
		  
		  feature="nspikesInBurst"; print ("Running nSpikes in bursts distribution analysis.")
		  params=parameters$burst.distribution.nSpikes
		  p<-IGM.plot.distributions(s[[i]],minVals=params$min.cases,xlimit=params$x.lim,binsInSec=params$bins.in.seg,feature=feature,
		                        filterValuesByMin=params$filter.by.min,minValues=params$min.values,perWell=params$per.well,outputdir=outputdir)}
		
		# Perform duration of bursts distribution analysis
		if (parameters$burst.distribution.durn$perform){
		  
		  feature="duration"; print ("Running duration of bursts distribution analysis.")
		  params=parameters$burst.distribution.durn
		  p<-IGM.plot.distributions(s[[i]],minVals=params$min.cases,xlimit=params$x.lim,binsInSec=params$bins.in.seg,feature=feature,
		                        filterValuesByMin=params$filter.by.min,minValues=params$min.values,perWell=params$per.well,outputdir=outputdir)}
		
		
		# Perform duration of bursts distribution analysis
		if (parameters$burst.distribution.spikeFreq$perform){
		  
		  feature="spikesDensityInBurst"; print ("Running spike density in bursts distribution analysis.")
		  params=parameters$burst.distribution.spikeFreq
		  p<-IGM.plot.distributions(s[[i]],minVals=params$min.cases,xlimit=params$x.lim,binsInSec=params$bins.in.seg,feature=feature,
		                        filterValuesByMin=params$filter.by.min,minValues=params$min.values,perWell=params$per.well,outputdir=outputdir)}
		
		#MFR
		p<- .plot.meanfiringrate(s[[i]], main = "Mean Firing Rate by Plate (Hz)")
		#p<- plot(s[[i]], main = "", label.cells = FALSE, use.names = FALSE)

  		p<-.channel.plot.by.well(s[[i]],resp="meanfiringrate", resp.label="Mean Firing Rate (Hz)")
		#Mean Duration
  		p<-.channel.plot.by.well(s[[i]],resp="bs$mean.dur",resp.label="Mean Duration of Burst (s)")
		#plot of Number of bursts by channel and well  
  		p<-.channel.plot.by.well(s[[i]],resp="bs$nbursts",resp.label="Number of Bursts")
		#mean Inter Burst Interval 
  		p<-.channel.plot.by.well(s[[i]],resp="bs$mean.IBIs",resp.label="Mean IBIs (ms)")
		# mean ISI within bursts 
  		p<-.channel.plot.by.well(s[[i]],resp="bs$mean.isis",resp.label="Mean ISI w/i Bursts (s)")
		#mean burst per minute 
  		p<-.channel.plot.by.well(s[[i]],resp="bs$bursts.per.min",resp.label="Mean Burst per Minute")
		#mean spikes in a burst
  		p<-.channel.plot.by.well(s[[i]],resp="bs$mean.spikes",resp.label="Mean # Spikes/Burst")
		#% spikes in a burst
  		p<-.channel.plot.by.well(s[[i]],resp="bs$per.spikes.in.burst",resp.label="% Spikes/Burst")
		dev.off()  
	}
}

write.plate.summary.for.bursts<-function(s,outputdir) {
	masterSum<-.get.burst.info.averaged.over.well(s)
	csvwell <- paste(outputdir,"/",get.project.plate.name(s[[1]]$file),"_well_bursts.csv",sep="")

	for (i in 1:length(s)){
		div <- .get.div(s[[i]])
		basename <- get.file.basename(s[[i]]$file)
		csvfile<-paste(outputdir,"/",basename,"_bursts.csv",sep="")

	  	##########data frame summarized over well
  		#get number of object in masterSum[[1]] list
  		tempdf<-c(); tempcolnames<-c()
  		for (j in 2:length(masterSum[[i]])){
    			tempc<-unlist(masterSum[[i]][j])
    			tempdf<-cbind(tempdf,tempc)
    			tempcolnames<-c(tempcolnames,names(masterSum[[i]][j]) )
    		}#end of loop through masterSum list objects
  
  		#need to switch around columns so first columns come first
		if (dim(tempdf)[2] > 20) { #for now
			if (dim(tempdf)[1] == 1) {
				df<-cbind(t(tempdf[,21:25]),t(tempdf[,1:20]))
			} else {
    				df<-cbind(tempdf[,21:25],tempdf[,1:20])
			}
    			colnames<-c(tempcolnames[21:25],tempcolnames[1:20])
    			colnames(df)<-colnames
		}

  		##################channel by channel burst summary

  		#meta data and misc
  		#get vector of which wells are active
  		wellindex<-which(is.element(names(s[[i]]$treatment), unique(s[[i]]$cw)) )
  		well<-c(); treatment<-c(); size<-c(); dose<-c(); file<-c();
  
  		file<-rep(strsplit( basename(s[[i]]$file),".RData")[[1]][1], length(s[[i]]$cw))
  		well<-s[[i]]$cw
  
  		#channel data frame
  		df2<-cbind(file,well,as.data.frame( s[[i]]$bs[1:length(s[[i]]$bs)]) )
  
 		
      	#write a title
    		write.table("Bursting Analysis averaged over Each Well",
                csvfile, sep=",", append=FALSE,row.names=FALSE,col.names=FALSE) 
  

      	write.table(paste("file= ",  strsplit( basename(s[[i]]$file),".RData")[[1]][1], sep=""),
                  csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE) 
      	write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
      	#recording time
      	write.table(paste("recording time (s): [", paste(s[[i]]$rec.time[1],round(s[[i]]$rec.time[2]), sep=" ,"),
            "]",sep=""),csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
  
      	write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
  
      	#summary write data, Sahar 26/11/2014 - add back genotype (column 1)
		if (dim(df)[1] == 1) {
			suppressWarnings(write.table(t(df[,-c(2:3)]),
                  csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=TRUE))
			suppressWarnings(write.table(cbind(div,t(df[,-c(2:3)])),
                  csvwell, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE))
		} else {
			suppressWarnings(write.table(df[,-c(2:3)],
                  csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=TRUE))
			suppressWarnings(write.table(cbind(div,df[,-c(2:3)]),
                  csvwell, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE))
		}
      	
      	#new lines
      	write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
      	write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
      
      	#title
      	write.table("Channel Burst Summary",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
      	write.table(paste("file= ",  strsplit( basename(s[[i]]$file),".RData")[[1]][1], sep=""),
                  csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE) 
      	write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
      
      	#channel data
      	suppressWarnings(write.table(df2,
                  csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=TRUE))
    		#new lines
    		write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
    		write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
  		write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
  		write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
  		write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
  		write.table(" ",csvfile, sep=",", append=TRUE,row.names=FALSE,col.names=FALSE)
	}#end of loop through writting tables
}
