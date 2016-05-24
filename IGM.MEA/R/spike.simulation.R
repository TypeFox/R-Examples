##read the data and have access to all the meta-data
.spike.simulation<-function(s1,
                           elec.min.rate=(1/60), 
                           elec.max.rate=25,
                           well.min.rate=15){
  dt <- 1
  current.electrode.rate <- s1$meanfiringrate
  rec.start <- s1$rec.time[1]
  rec.end <- s1$rec.time[2]
  spikes <- list()
  for (electrode in 1:length(s1$spikes)) {
    rate <- current.electrode.rate[electrode] * dt / 1000.0
    spiketimes <- list()
    timepoints <- seq(rec.start,rec.end, by = 0.001)
    p <- which(rate > runif(length(timepoints)))
    spikes[[electrode]] <- timepoints[which(rate > runif(length(timepoints)))]
    
  }
  names(spikes) <- names(s1$spikes)
  temp.s <- construct.s(spikes, NULL, time.interval = 1, beg = NULL, end = NULL, corr.breaks = 0, 
                        s1$layout, filename = s1$file)
  
  
  #indices of low and high firing rate
  low <- which(temp.s$meanfiringrate < elec.min.rate)
  high <- which(temp.s$meanfiringrate > elec.max.rate)
  
  ## TODO, check that low and high are non-zero length vectors.
  extremes <- c(low, high)
  
  bad.ids <- names(extremes)
  bad.ids <- c("-", bad.ids)  # "-" needed to remove these ids!
  
  s2 <- remove.spikes(temp.s, bad.ids)
  
  s2$treatment<-s1$treatment
  s2$size<-s1$size
  s2$units<-s1$units
  s2$dose<-s1$dose
  s2$well<-s1$well
  
  #get.num.AE
  s2<-get.num.AE(s2)
  
  #indices of low and high firing rate
  
  low <- which(s2$nAE < well.min.rate)
  
  bad.wells <- names(low)
  bad.wells <- c("-", bad.wells)   # "-" needed to remove these well!
  #just these three for example
  s<- remove.spikes(s2, bad.wells)
  
  s$goodwells<-names(which(s2$nAE >= well.min.rate))
  
  #[which(s2$nAE >= well.min.rate)
  s$treatment<-s1$treatment
  names(s$treatment)<-s1$well
  s$size<-s1$size
  names(s$size)<-s1$well
  s$units<-s1$units
  names(s$units)<-s1$well
  s$dose<-s1$dose
  names(s$dose)<-s1$well
  s$well<-s1$well
  s<-get.num.AE(s)
  s$timepoint <- s1$timepoint
  if (s$nspikes[1] >0) {
    s$allb <- lapply(s$spikes, mi.find.bursts,s$parameters$mi.par)
    s$bs<-calc.burst.summary(s)
  }
  
  s <- calculate.isis(s)
  s$well.stats <- IGM.mean.firingrate.by.well(s)
  s
}
