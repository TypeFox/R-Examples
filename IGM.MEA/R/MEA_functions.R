
get.data<-function(caption=""){

  #get the directory containing the .spikelist files
  spikeFiles<-sort(tk_choose.files(caption=caption) )
   
  return(spikeFiles)  
}

get.num.AE<-function(s2){
  	#add number of active electrodes
  	s2$nAE<-rep(0,length(s2$well))
  	names(s2$nAE)<-s2$well
  	for (i in 1:s2$NCells){
    		s2$nAE[which(substr(s2$channels[i],1,2)==(s2$well))]=
      	s2$nAE[which(substr(s2$channels[i],1,2)==(s2$well))]+1
    
    		s2$cw[i]<-substr(s2$channels[i],1,2)
  	}
  	s2 
}

#this function removes bad channels
remove.spikes <- function(s, ids) {
  	## ids=vector of indicies eg c(1,2,4,5)
  	## Remove spikes listed in IDS from S data structure, and return
  	## new structure.
  
  	beg <- s$rec.time[1]
  	end <- s$rec.time[2]
  	corr.breaks <- 0         #TODO: hardcoded for axion!
  	layout <- s$layout
  	filename <- s$file#paste0(s$file, ".edited")
  	s2 <- construct.s(s$spikes, ids, s$rates$time.interval, beg, end,
                    corr.breaks, layout, filename)
  	s2
}

################THIS FUNCTION NEEDS EDITING TO PASS BACK WELL INFO#####
#this function runs through spikes, cutting out
#bad wells and bad electrodes, it also adds the well data

#this function returns a list with the chemical names for corresponding
#date, plate number and wells found in "file"
chem.info.2<-function(file,masterChemFile=masterChemFile) {
  	#get chemical info
  	#*********load chemical list*************
  	# read in masterChemFile
  	masterChem=read.csv(masterChemFile)
  	masterCD<-as.data.frame(masterChem)
  	#remove extraneous NA columns
  	masterCD<-masterCD[,1:9]
  
  	#reconstruct file names
  	temp1<-paste(masterCD$Project, masterCD$Experiment.Date,masterCD$Plate.SN,sep="_")
  
  	#add a new column called masterCD$filenames
  	masterCD$filenames<-temp1
  
  	#ensure log file is ordered correctly so it can be read in correctly
  	masterCD<-masterCD[order(masterCD$Experiment.Date,masterCD$Plate.SN,masterCD$Well),]
  
  	#****match wells to chemicals *****
  	shortFileName<-paste( strsplit(basename(file),"_")[[1]][1],
                        strsplit(basename(file),"_")[[1]][2],
                        strsplit(basename(file),"_")[[1]][3],sep="_")
  
  	plate.chem.info<-list()
  	count=1;
    matchedFileName = 0;
  	for (i in which(shortFileName==masterCD$filename) ){
  	  matchedFileName = 1;
    		#get all info from chem list
    		plate.chem.info$well[count]<-paste(masterCD$Well[i])
    		plate.chem.info$treatment[count]<-paste(masterCD$Treatment[i])
    		plate.chem.info$size[count]<-paste(masterCD$Size[i])
    		plate.chem.info$dose[count]<-paste(masterCD$Dose[i])
    		plate.chem.info$units[count]<-paste(masterCD$Units[i])
    		count=count+1
    
  	}#end of for loop through masterCD$file
    if (matchedFileName == 0){
      print(paste("File ",shortFileName," was not found in the possible file names 
                  constructed from exp log file:",unique(masterCD$filename),sep=""))
    }
  	if (!is.element(length(plate.chem.info$well),c(12,48)) ){
    		print(paste("Info exists for ",length(plate.chem.info$well),
                " wells; Some wells have no data.", sep=""))
  	}
  	plate.chem.info
}

#purpose: given a list containing spikes and s$bs containing burst info,
# plots the resp=response variable, by channel, in a lattice grid grouped by well
#output=a plot handle, p
#input: spikes and respsonse variable
#EXAMPLE:    p<-.channel.plot.by.well(s,resp="meanfiringrate")
#EXAMPE: list nested in list: p<-.channel.plot.by.well(s,resp="bs$mean.dur")
.channel.plot.by.well<-function(s , resp, resp.label ){   
  	par(mfrow=c(1,1))  
  	if (length(s$well)<=12){
    		well.layout=c(4,3)
    		well.names <- paste(rep(LETTERS[3:1], each = 4), rep(1:4, 3), sep = "")
    		treatment_size<-paste(c(s$treatment[9:12],s$treatment[5:8],s$treatment[1:4]),
                          c(s$size[9:12],s$size[5:8],s$size[1:4]),sep=" ")
    		names(well.names) <- paste( paste(rep(LETTERS[3:1], each = 4), rep(1:4, 3), sep = ""),
                                treatment_size,sep='=')
    		par.strip = list(cex = 1) 

  	} else {
    		well.layout=c(8,6)   
    		well.names<-paste(rep(LETTERS[6:1], each = 8), rep(1:8, 6), sep = "")
    		treatment_size<-paste(c(s$treatment[41:48],s$treatment[33:40],s$treatment[25:32],
                            s$treatment[17:24],s$treatment[9:16],s$treatment[1:8]),
                          c(s$size[41:48],s$size[33:40],s$size[25:32],
                            s$size[17:24],s$size[9:16],s$size[1:8]),sep=" ")
    		names(well.names) <- paste( paste(rep(LETTERS[3:1], each = 4), rep(1:4, 3), sep = ""),
                                treatment_size,sep='=')
    		par.strip = list(cex = .6) 
    
  	}
  
  	s$active.wells<-axion.elec2well(s$channels)
  
  
  	if (length(strsplit(resp,"$",fixed=TRUE)[[1]])>1 ){
    		response<-get(strsplit(resp,"$",fixed=TRUE)[[1]][2] , get(strsplit(resp,"$",fixed=TRUE)[[1]][1], s) )
  	} else {
    		response<-get(strsplit(resp,"$",fixed=TRUE)[[1]][1],  s)
  	}
  
  	p <- xyplot(response ~ factor(channels) | 
                factor(active.wells, labels=names(well.names),levels = well.names),
              data = s, drop.unused.levels = FALSE, layout = well.layout, 
              xlab = "Channels within well",
              ylab = paste(resp.label,sep=""), pch = 20 ,
              main= paste( paste(resp.label, " by Channels within Wells",sep=""), 
                                  paste("file= ",  strsplit( basename(s$file),".RData")[[1]][1], sep=""),                             
                                  sep='\n') , 
              scales = list(x = list(relation = "free",
                                     draw = FALSE)),
              par.strip.text = par.strip )

  
  	print( p)
  	p
}#end of .channel.plot.by.well  

#********average and sum burst variables across each well
#input: s is a list containing burst info, and meta data
##purpose: average across wells
##the list returned, (masterSum[[1]],masterSum[[2]]..etc for each spike list s[[i]])
##has meta data and has been filtered according to weather it's 48 or 12 well
#necessary: the timepoint "00" is needed to set which wells are active etc
.get.burst.info.averaged.over.well<-function(s){  
  	masterSum<-list()#summary over all files
  	for (i in 1:length(s)){
    		sum=list()#summary for each timepoint
    		#calculate bursting variables for current data File
    		nbursts <- sapply(s[[i]]$allb, nrow) 
    		allb<-s[[i]]$allb
    		tempsum <- calc.burst.summary(s[[i]])
    		
    		#ISIs: gets the ISI for each channel of s[[i]]
    		ISIs = .calc.all.isi(s[[i]], allb)
    
    		#IBIs get IBI's across all inter burst intervals across all data
    		tempIBIs<-.calc.all.ibi(s[[i]],allb)
    
    		#loop through goodwells
    		for (j in 1:length(s[[i]]$goodwells)){
      		#indicator of current well
        		icurrentwell<-(s[[i]]$goodwells[j]==s[[i]]$cw)
        
        		#index of current well
        		incurrentwell<-which(s[[i]]$goodwells[j]==s[[i]]$cw)
        
        		if (sum(icurrentwell)!=0){       
        			#####variables that need summing and averaging  
        			#total spikes across all AE in current well
        			sum$nspikes[j]<-sum(tempsum$spikes[icurrentwell], na.rm=TRUE)
        			sum$nAB[j] <- length(which(nbursts[incurrentwell]>0))
        			#Total recorded time on current well= recording time * nAE
        			sum$duration[j]<-length(incurrentwell)*(s[[i]]$rec.time[2]-s[[i]]$rec.time[1])
        
        			#mean duration
        			sum$mean.dur[j]<-mean(tempsum$mean.dur[incurrentwell],na.rm=TRUE)
        
        			#mean spikes per second
        			sum$mean.freq[j]<-mean(tempsum$mean.freq[incurrentwell],na.rm=TRUE)
        			#total number of bursts
        			sum$nbursts[j]<-sum(tempsum$nbursts[icurrentwell], na.rm=TRUE)
        			#mean burst per second
        			sum$bursts.per.sec[j]<-mean(tempsum$bursts.per.sec[incurrentwell])
        			#mean burst per minute
        			sum$bursts.per.min[j]<-sum$bursts.per.sec[j]*60
        
        			#finds the mean of duration for a particular well (across all channels)
        			#burstinfo(allb[icurrentwell],"durn") takes out the column "durn" of all
        			#matricies allb among the indicator set icurrentwell
        			#get duration data across all channels of current well
        			sum$mean.dur[j]<-mean(unlist(burstinfo(allb[icurrentwell],"durn")),na.rm=TRUE)
        
        			#sd of current well burst durations
        			sum$sd.dur[j]<-sd(unlist(burstinfo(allb[icurrentwell],"durn")))
        
        			#mean frequency within a burst
        			sum$mean.freq.in.burst[j]<-
          				mean(unlist(burstinfo(allb[incurrentwell],"len"))/
                 			unlist(burstinfo(allb[incurrentwell],"durn")), na.rm=TRUE)
        
        			#sd frequency within a burst
        			sum$sd.freq.in.burst[j]<-sd(unlist(burstinfo(allb[incurrentwell],"len"))/
                                      unlist(burstinfo(allb[incurrentwell],"durn")),na.rm=TRUE)
        
        			#mean of ISI across all channels in current well
        			sum$mean.ISIs[j] = mean(unlist(ISIs[incurrentwell]), na.rm=TRUE)
        
        			#finds sd of ISI across all channels in current well
        			sum$sd.ISIs[j] = sd( unlist(ISIs[incurrentwell]), na.rm = TRUE)
        
        			#len=#spikes in burst (length of burst in bursts)
        			#mean.spikes.in.burst
        			ns<-unlist(burstinfo(allb[icurrentwell],"len"))   
        			sum$mean.spikes.in.burst[j] <- round(mean(ns,na.rm=TRUE), 3)
        
        			#sd of spikes in burst
        			sum$sd.spikes.in.burst[j] <- round(sd(ns,na.rm=TRUE), 3)
        
        			#total number of spikes arcross all bursts
        			sum$total.spikes.in.burst[j] <- sum(ns, na.rm=TRUE)
        
        			#percent of spikes in bursts
        			sum$per.spikes.in.burst[j] <- 
          				round(100 * (sum$total.spikes.in.burst[j]/sum$nspikes[j]),3)
        
        			#mean IBI
        			sum$mean.IBIs[j]<-round(mean(unlist(tempIBIs[incurrentwell]),na.rm=TRUE),3)
        			#sd IBI
        			sum$sd.IBIs[j]<-round(sd(unlist(tempIBIs[incurrentwell]),na.rm=TRUE),3)
        			#cv IBI
        			sum$cv.IBIs[j]<-round(sum$mean.IBIs[j]/sum$sd.IBIs[j],3)
				
      		} else {
        			sum$nspikes[j]<-NA
				sum$nAB[j] <- NA
        			sum$duration[j]<-NA
        			sum$mean.dur[j]<-NA
        			sum$mean.freq[j]<-NA
        			sum$nbursts[j]<-NA
        			sum$bursts.per.sec[j]<-NA
        			sum$bursts.per.min[j]<-NA
                		sum$mean.dur[j]<-NA
        			sum$sd.dur[j]<-NA
        			sum$mean.freq.in.burst[j]<-NA
        			sum$sd.freq.in.burst[j]<-NA
        			sum$mean.ISIs[j] = NA
        			sum$sd.ISIs[j] = NA
        			sum$mean.spikes.in.burst[j] <- NA
        			sum$sd.spikes.in.burst[j] <- NA
        			sum$total.spikes.in.burst[j] <- NA
        			sum$per.spikes.in.burst[j] <- NA
        			sum$mean.IBIs[j]<-NA
        			sum$sd.IBIs[j]<-NA
        			sum$cv.IBIs[j]<-NA

      		}# end of if/else for channels that dropped off
      
    		}#end of loop through goodwells
    
    		###Set all names
    		for (k in 1:length(names(sum))){
      		names(sum[[k]])=s[[i]]$goodwells
    		} 
    
    		#make a masterSum, that is a list of all the summaries
    		goodwellindex<-which(is.element(s[[i]]$well, s[[i]]$goodwells ) )
    
    		masterSum[[i]]<-sum
    		masterSum[[i]]$file<-strsplit(basename(s[[i]]$file),".RData")[[1]][1]
    		masterSum[[i]]$treatment<-s[[i]]$treatment[goodwellindex]
    		masterSum[[i]]$size=s[[i]]$size[goodwellindex]
    		masterSum[[i]]$dose=s[[i]]$dose[goodwellindex]
    		masterSum[[i]]$well<-s[[i]]$well[goodwellindex]
    		masterSum[[i]]$nAE<-s[[i]]$nAE[goodwellindex]
    		masterSum[[i]]$timepoint=rep(s[[i]]$timepoint[1],length(s[[i]]$goodwells))
    		masterSum[[i]]$start.rec.time<-rep( s[[i]]$rec.time[1],length(s[[i]]$goodwells) )
    		masterSum[[i]]$end.rec.time<-rep( s[[i]]$rec.time[2],length(s[[i]]$goodwells) )
    		masterSum[[i]]$goodwells<-s[[i]]$goodwells
    
  	}#end of for loop through sum/averaging burst variables
  
  	masterSum
}#end of get.burst.info
