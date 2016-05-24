# burst_functions.R

#this function takes in a list of spikes (s), with burst info
# digitizes spikes into 10ms bins (1 in bin if any spike, 0 otherwise)
# then computes correlation btw ae above, below, to left and to right,
# and at 4 corners. Then finds the mean corr. for each electode
# input: s[[i]] object
# returns: plate.corr: a mean correlation for each active channel
#  between the neighboring electrodes.
#  if no neurons are neighboring, NA is returned in channel
.local.corr.by.channel<-function(s) {
  
  # store data
  plate.corr=list(plate.name=list(), local.cor=list())
  
  #loop through different time points
  #t.p is time points
  for (t.p in c(1:length(s))){
    
    #now get correlations within well
    temp=list()
    r.vector=list()
    r.vector[[length( unique(s[[t.p]]$cw) )]]=0
    
    names(r.vector)=unique(s[[t.p]]$cw)
    
    #using  (Demas 2009), figure 5, closer neurons are more correlated
    #therefore I'm going to measure only cor. btw adjacent electrodes
    
    #loop through wells
    for (w in c(1:length(unique(s[[t.p]]$cw))) ){
      
      #get indices of current well
      index.w=which(s[[t.p]]$cw==unique(s[[t.p]]$cw)[w])
      
      #convert seconds to ms
      start.ms=round(s[[t.p]]$rec.time[1]*10^2,digits=0)
      stop.ms=round(s[[t.p]]$rec.time[2]*10^2,digits=0)
      total.ms=stop.ms-start.ms
      
      
      #digitize spikes
      dig.s=list()
      dig.s[[length(index.w)]]=rep(0, total.ms )
      
      
      #digitize spikes: loop through channels
      for (i in c(1:length(index.w))){
        
        dig.s[[i]]=rep(0, total.ms)
        current.s.ms=round(s[[t.p]]$spikes[[index.w[i] ]]*10^2, digits=0)- start.ms 
        for (j in current.s.ms ){
          #get millisecond
          dig.s[[i]][ j ] = 1
        }
      }
      
      #get indices of current well
      index.c=which(s[[t.p]]$cw==unique(s[[t.p]]$cw)[w])
      
      #get channels
      w.rows=as.numeric( substr(s[[t.p]]$channels[index.c],4,4) )
      w.cols=as.numeric( substr(s[[t.p]]$channels[index.c],5,5) )
      coords=rbind(w.rows,w.cols)
      
      #initialize c.vector: vector of mean(r), for each channel in current well
      c.vector = c()
      
      #loop through channels, getting cor between an electrode
      #and electrodes above, below, to left and to right
      # if those electrodes have AE
      for (i in c(1:length(s[[t.p]]$channels[index.c] ) ) ){
        
        #c.c is current coordinates of electrode on plate map
        c.c=coords[,i]
        temp.corr= c()
        
        #check for suitable coordinates to coorelate with       
        #example: take c.c = (4,4)
        # (c.c[1]==coords[1,] & c.c[2]+1==coords[2,]) means:
        # is current channel row (c.c[1]=4) same as and row in other electodes,
        # and current column (c.c[2]=4), +1, =5 the same as any column
        # this would correspond to an electrode (4,5). none exists so if loop not entered
        
        if (sum(c.c[1]==coords[1,] & c.c[2]+1==coords[2,])>0){
          # (4,5)
          ind.c.to.comp=which(c.c[1]==coords[1,] & c.c[2]+1==coords[2,])
          temp.corr=c(temp.corr, cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum( c.c[1]==coords[1,] & c.c[2]-1==coords[2,])>0){
          #this would be (4,3),
          ind.c.to.comp= which(c.c[1]==coords[1,] & c.c[2]-1==coords[2,]) 
          temp.corr=c(temp.corr, cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum(c.c[1]-1==coords[1,] & c.c[2]==coords[2,])>0){
          # this would be (3,4), so no match in 48 well plate
          ind.c.to.comp= which(c.c[1]-1==coords[1,] & c.c[2]==coords[2,]) 
          temp.corr=c(temp.corr, cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum(c.c[1]+1==coords[1,] & c.c[2]==coords[2,])>0){
          # this would be (5,4)
          ind.c.to.comp=which(c.c[1]+1==coords[1,] & c.c[2]==coords[2,]) 
          temp.corr=c(temp.corr, cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum(c.c[1]+1==coords[1,] & c.c[2]+1==coords[2,])>0){
          # this would be (5,5)
          ind.c.to.comp=which(c.c[1]+1==coords[1,] & c.c[2]+1==coords[2,]) 
          temp.corr=c(temp.corr, cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum(c.c[1]-1==coords[1,] & c.c[2]-1==coords[2,])>0){
          # this would be (3,3)
          ind.c.to.comp=which(c.c[1]-1==coords[1,] & c.c[2]-1==coords[2,]) 
          temp.corr=c(temp.corr, cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum(c.c[1]+1==coords[1,] & c.c[2]-1==coords[2,])>0){
          # this would be (5,3)
          ind.c.to.comp=which(c.c[1]+1==coords[1,] & c.c[2]-1==coords[2,]) 
          temp.corr=c(temp.corr, cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum(c.c[1]-1==coords[1,] & c.c[2]+1==coords[2,])>0){
          # this would be (3,5)
          ind.c.to.comp=which(c.c[1]-1==coords[1,] & c.c[2]+1==coords[2,]) 
          temp.corr=c(temp.corr, cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        #vector of correlations inside well
        c.vector = c(c.vector, mean(temp.corr, na.rm=T ) )
        
      }#end of loop through channels
      names(c.vector)= s[[t.p]]$channels[index.c]
      r.vector[[w]] = c.vector
      
    }#end of loop through wells
    
    plate.corr$plate.name[t.p]=strsplit(basename(s[[t.p]]$file),'.RData')[[1]][1]
    
    plate.corr$local.cor[[t.p]]=r.vector
    
    
    
  }#end of loop through time points
  plate.corr # returns plate.corr
  
} # end of .local.corr










#this function takes in a list of spikes (s), with burst info
# digitizes spikes into 10ms bins (1 in bin if any spike, 0 otherwise)
# then computes correlation btw ae above, below, to left and to right
# input: s[[i]] object
# returns: plate.corr
.local.corr<-function(s) {
  
  # store data
  plate.corr=list(plate.name=list(), local.cor=list())
  
  #loop through different time points
  #t.p is time points
  for (t.p in c(1:length(s))){
    
    
    
    #now get correlations within well
    temp=list()
    r.vector=list()
    r.vector[[length( unique(s[[t.p]]$cw) )]]=0
    
    names(r.vector)=unique(s[[t.p]]$cw)
    
    #using  (Demas 2009), figure 5, closer neurons are more correlated
    #therefore I'm going to measure only cor. btw adjacent electrodes
    
    #loop through wells
    for (w in c(1:length(unique(s[[1]]$cw))) ){
      
      #get indices of current well
      index=which(s[[t.p]]$cw==unique(s[[1]]$cw)[w])
      
      #convert seconds to ms
      start.ms=round(s[[t.p]]$rec.time[1]*10^2,digits=0)
      stop.ms=round(s[[t.p]]$rec.time[2]*10^2,digits=0)
      total.ms=stop.ms-start.ms
      
      
      #digitize spikes
      dig.s=list()
      dig.s[[length(index)]]=rep(0, total.ms )
      
      
      #digitize spikes: loop through channels
      for (i in c(1:length(index))){
        
        dig.s[[i]]=rep(0, total.ms)
        current.s.ms=round(s[[t.p]]$spikes[[index[i] ]]*10^2, digits=0)- start.ms 
        for (j in current.s.ms ){
          #get millisecond
          dig.s[[i]][ j ] = 1
        }
      }
      
      
      
      #get indices of current well
      index=which(s[[t.p]]$cw==unique(s[[1]]$cw)[w])
      
      #get channels
      w.rows=as.numeric( substr(s[[t.p]]$channels[index],4,4) )
      w.cols=as.numeric( substr(s[[t.p]]$channels[index],5,5) )
      coords=rbind(w.rows,w.cols)
      
      #loop through channels, getting cor between an electrode
      #and electrodes above, below, to left and to right
      # if those electrodes have AE
      for (i in c(1:length(s[[1]]$channels[index]))){
        c.c=coords[,i]
        
        #check for suitable coordinates to coorelate with
        if (sum(c.c[1]==coords[1,] & c.c[2]+1==coords[2,])>0){
          ind.c.to.comp=which(c.c[1]==coords[1,] & c.c[2]+1==coords[2,])
          r.vector[[w]]=c(r.vector[[w]], cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum( c.c[1]+1==coords[1,] & c.c[2]-1==coords[2,])>0){
          ind.c.to.comp= which(c.c[1]+1==coords[1,] & c.c[2]-1==coords[2,]) 
          r.vector[[w]]=c(r.vector[[w]], cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum(c.c[1]+1==coords[1,] & c.c[2]==coords[2,])>0){
          ind.c.to.comp= which(c.c[1]+1==coords[1,] & c.c[2]==coords[2,]) 
          r.vector[[w]]=c(r.vector[[w]], cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        if (sum(c.c[1]+1==coords[1,] & c.c[2]+1==coords[2,])>0){
          ind.c.to.comp=which(c.c[1]+1==coords[1,] & c.c[2]+1==coords[2,]) 
          r.vector[[w]]=c(r.vector[[w]], cor(dig.s[[i]], dig.s[[ind.c.to.comp]]) )
        }
        
      }#end of loop through channels
      
    }#end of loop through wells
    
    plate.corr$plate.name[t.p]=strsplit(basename(s[[t.p]]$file),'.RData')[[1]][1]
    
    plate.corr$local.cor[[t.p]]=r.vector
    
    
    
  }#end of loop through time points
  plate.corr # returns plate.corr
  
} # end of .local.corr












#this function takes in a list of spikes (s), with burst info
# digitizes spikes into 10ms bins (1 in bin if any spike, 0 otherwise)
# then computes correlation btw that electrode and every other electrode.
# Then finds the mean corr. for each electode, and returns a mean r for each electrode
# input: s[[i]] object
# returns: plate.corr: a mean correlation for each active channel
#  between the neighboring electrodes.
#  if no neurons are neighboring, NA is returned in channel
.local.corr.all<-function(s) {
  
  # store data
  plate.corr=list(plate.name=list(), local.cor=list())
  
  #loop through different time points
  #t.p is time points
  for (t.p in c(1:length(s))){
    
    
    
    temp=list()
    #r.vector is a list for each channel
    r.vector=list()
    r.vector[[length( unique(s[[t.p]]$cw) )]]=0
    
    names(r.vector)=unique(s[[t.p]]$cw)
    
    #using  (Demas 2009), figure 5, closer neurons are more correlated
    #therefore I'm going to measure only cor. btw adjacent electrodes
    
    #loop through wells
    for (w in c(1:length(unique(s[[t.p]]$cw))) ){
      
      #get indices of channels within current well
      index.w=which(s[[t.p]]$cw==unique(s[[t.p]]$cw)[w])
      
      #convert seconds to ms
      start.ms=round(s[[t.p]]$rec.time[1]*10^2,digits=0)
      stop.ms=round(s[[t.p]]$rec.time[2]*10^2,digits=0)
      total.ms=stop.ms-start.ms
      
      
      #digitize spikes
      dig.s=list()
      dig.s[[length(index.w)]]=rep(0, total.ms )
      
      
      #digitize spikes: loop through current channels
      #dig.s[[i]] is a list 1:# channels within current well
      for (i in c(1:length(index.w))){
        
        dig.s[[i]]=rep(0, total.ms)
        current.s.ms=round(s[[t.p]]$spikes[[index.w[i] ]]*10^2, digits=0)- start.ms 
        for (j in current.s.ms ){
          #get millisecond
          dig.s[[i]][ j ] = 1
        }
      } #end loop through current channels 
      
      
      #get indices of current well
      index.c=which(s[[t.p]]$cw==unique(s[[t.p]]$cw)[w])
      
      #initialize c.vector: vector of mean(r), for each channel in current well
      c.vector = c()
      
      #loop through channels in current well, getting cor between an electrode
      #and electrodes above, below, to left and to right
      # if those electrodes have AE
      for (curCh in c(1:length(s[[t.p]]$channels[index.c] ) ) ){
        
        #c.c is current coordinates of electrode on plate map
        temp.corr= c()
        
        #loop through channels other than current one
        for (otherCh in setdiff(c(1:length(s[[t.p]]$channels[index.c] )),curCh)  ){
          temp.corr=c(temp.corr, cor(dig.s[[curCh]], dig.s[[otherCh]]) )
          
        } #end of loop through other channels
        
        
        #vector of mean correlations inside well by channel
        c.vector = c(c.vector, mean(temp.corr, na.rm=T ) )
        
      }#end of loop through channels
      names(c.vector)= s[[t.p]]$channels[index.c]
      # each well get one vector
      r.vector[[w]] = c.vector
      
    }#end of loop through wells
    
    plate.corr$plate.name[t.p]=strsplit(basename(s[[t.p]]$file),'.RData')[[1]][1]
    
    plate.corr$local.cor[[t.p]]=r.vector
    
    
    
  }#end of loop through time points
  plate.corr # returns plate.corr
  
} # end of .local.corr.all





# get network spikes. 
# input: an s object
# output: 1. ns$measures[,1] give times of all ns 
#   2. ns$measures[,3] gives peak.val aka n spikes in ns for all ns
#    ns$measures[,4] gives durn=durn of ns for each ns
#   ns$brief$n gives total number of ns

.plate.ns<-function(s, t.p){
  

    
    ns=list()
    ns[[ length( s[[t.p]]$goodwells ) ]] = 0
    for (well in c(1:length( s[[t.p]]$goodwells ) ) ){
      
      ns[[well]] <- .compute.ns(s[[t.p]], ns.T=0.003, ns.N=6, 
                               sur=100, whichcells=s[[t.p]]$goodwells[well]  )
      ns[[well]]$well = paste( s[[t.p]]$goodwells[well] )
      print (paste("well ", paste(s[[t.p]]$goodwells[well]), " finished" )  )
      
    } #end of for loop through wells in a plate
    ns
  
} #end of .plate.ns


# this function calculates the rate of a timestamp variable at each minute, computes sd
# then gives the cv_time and cv_network
# input: 1) timestamp= list(timestamp[[channel1]],...,timestamp[[channelN]])
#         3) s = the objects
#  eg inputs:  s
#              timestamps = burst.train 
.cv.timestamps<-function(timestamps, s, t.p){
  
  #for (t.p in c(1:length(s)) ){

  
  
  #loop through wells
  cv_time=c() #cv_time for each well
  cv_network = c()  #cv_network for each network
  for (w in c(1:length(unique(s[[t.p]]$cw) ) ) ){
    
    #get indices of current well
    index=which(s[[t.p]]$cw==unique(s[[t.p]]$cw)[w])
    
    #convert seconds to minutes (m), and round up 
    total.m = round( (s[[t.p]]$rec.time[2] - s[[t.p]]$rec.time[1])/60, 0)  #total minutes
    
    #cv_time:
    #loop through channels in current well
    # cv_Ch is a vector for each well
    cv_ch = c() 
    for ( cur.ch in index ){
      #loop through minutes in current channel
      temp.mean=c()
      cur.beg = s[[t.p]]$rec.time[1]
      for (cur.min in c(1:(total.m-1) ) ){
        cur.end = cur.beg + 60 #look at 60 second window of time
        
        #n.min is used to get mean
        n.min = 1
        if ( cur.min == total.m-1 ){
          n.min = (s[[t.p]]$rec.time[2]%%60)/60 #this is modulus operator so remainder
        }
        #beg=first spike (out of total spikes) that's in burst
        #end=last spike (of nspikes) that's in burst     
        #IBI len = time between current burst and previous burst    
        #durn= duration of a burst  
        #mean.isis = interval between spikes in burst
        #SI = .surprise index
        # s[[1]]$spikes$F8_44[firstSpikeOfBurst], you can get time at which burst started by taking 'beg'
        # of current burst, firstSpikeOfBurst = s[[1]]$allb$F8_44[currentBurst,1]
        #per.spikes.out.burst mean.si   mean.isis sd.mean.isis mean.IBIs  sd.IBIs cv.IBIs
        temp.mean=c( temp.mean, 
                     sum( timestamps[[cur.ch]] <= cur.end & timestamps[[cur.ch]]>=cur.beg)  )/n.min
        cur.beg = cur.end #shift beginning of time window to end of last window
      } #end loop through minutes
      
      cv_ch = c( cv_ch, sd(temp.mean, na.rm=T )/mean(temp.mean, na.rm=T) )
      
    } # end of loop through index
    
    #cv_time is average cv well
    cv_time =  c(cv_time, mean( cv_ch, na.rm=T )  )
    
    
    
    
    
    
    #cv_network:
    #loop through channels in current well
    # cv_min is a vector for each well
    cv_min = c() 
    for (cur.min in c(1:(total.m-1) ) ) {
      #loop through minutes in current channel
      # temp.mean=c()
      cur.beg = s[[t.p]]$rec.time[1]
      cur.end = cur.beg + 60 #look at 60 second window of time
      
      
      #n.min is used to get mean
      n.min = 1
      if ( cur.min == total.m-1 ){
        n.min = (s[[t.p]]$rec.time[2]%%60)/60 #this is modulus operator so remainder
      }
      
      for ( cur.ch in index ){
        
        temp.mean=c( temp.mean, 
                     sum( timestamps[[cur.ch]] <= cur.end & timestamps[[cur.ch]]>=cur.beg)  )/n.min
        
      } #end loop through index
      
      cv_min = c( cv_min, sd(temp.mean, na.rm=T )/mean(temp.mean, na.rm=T) )
      
      cur.beg = cur.end #shift beginning of time window to end of last window
      
      
    } # end of loop through minutes
    
    #cv_time is average cv well
    cv_network =  c( cv_network, mean( cv_min, na.rm=T ) ) 
    
    
  }#end of loop through wells
  
  results = list( cv_time, cv_network)
  results
  
  # } #end of loop through time points
  
  
}  # end of .cv.timestamps function







# this function calculates the rate of a timestamp variable at each minute, computes sd
# then gives the cv_time and cv_network
# input: 1) timestamp= list(timestamp[[channel1]],...,timestamp[[channelN]])
#         3) s = the objects
#  eg inputs:  s
#              timestamps = burst.train 
.cv.feature<-function(feature, s){
  
  #for (t.p in c(1:length(s)) ){
  t.p=1
  
  
  #loop through wells
  cv_time=c() #cv_time for each well
  
  # loop through wells
  for (w in c(1:length(unique(s[[t.p]]$cw) ) ) ){
    
    #get indices of current well
    index=which(s[[t.p]]$cw==unique(s[[t.p]]$cw)[w])
    
    #convert seconds to minutes (m), and round up 
    total.m = round( (s[[t.p]]$rec.time[2] - s[[t.p]]$rec.time[1])/60, 0)  #total minutes
    
    #cv_time:
    #loop through channels in current well
    # cv_Ch is a vector for each well
    cv_ch = c() 
    for ( cur.ch in index ){
      
      #only go through the channel if there's burst data
      if (  !(  is.null( dim( feature[[cur.ch]] ) ) )  ){
        
        #loop through minutes in current channel
        temp.mean=c()
        cur.beg = s[[t.p]]$rec.time[1]
        for (cur.min in c(1:(total.m-1) ) ){
          cur.end = cur.beg + 60 #look at 60 second window of time
          
          
          #beg=first spike (out of total spikes) that's in burst
          #end=last spike (of nspikes) that's in burst     
          #IBI len = time between current burst and previous burst    
          #durn= duration of a burst  
          #mean.isis = interval between spikes in burst
          #SI = .surprise index
          # s[[1]]$spikes$F8_44[firstSpikeOfBurst], you can get time at which burst started by taking 'beg'
          # of current burst, firstSpikeOfBurst = s[[1]]$allb$F8_44[currentBurst,1]
          #per.spikes.out.burst mean.si   mean.isis sd.mean.isis mean.IBIs  sd.IBIs cv.IBIs
          index.t = (  feature[[cur.ch]][,1]<= cur.end & feature[[cur.ch]][,1]>=cur.beg)
          if (  sum(index.t)>0   ){
            
            temp.mean=c( temp.mean,  mean( feature[[cur.ch]][index.t, 2], na.rm=T)  ) 
          } # end if ( sum(index.t)>0  )
          cur.beg = cur.end #shift beginning of time window to end of last window
        } #end loop through minutes
        if (length(temp.mean)>1 ){
          cv_ch = c( cv_ch, sd(temp.mean, na.rm=T )/mean(temp.mean, na.rm=T) )
        } else{
          cv_ch = c( cv_ch, NA )
        }
        
      } else {
        
        cv_ch = c( cv_ch, NA )
      } #end of if (  !is.na(feature[[cur.ch]])  )
      
    } # end of cur.ch through index
    
    #cv_time is average cv well
    cv_time =  c(cv_time, mean( cv_ch, na.rm=T )  )
    
  }#end of loop through wells
  
  
  
  
  
  
  
  
  
  #cv_network:
  #loop through channels in current well
  # cv_min is a vector for each well
  cv_network = c()  #cv_network for each network
  for (w in c(1:length(unique(s[[t.p]]$cw) ) ) ){
    
    
    cv_min = c() 
    for (cur.min in c(1:(total.m-1) ) ) {
      #loop through minutes in current channel
      # temp.mean=c()
      cur.beg = s[[t.p]]$rec.time[1]
      cur.end = cur.beg + 60 #look at 60 second window of time
      
      
      for ( cur.ch in index ){
        
        
        
        if (  !(  is.null( dim( feature[[cur.ch]] ) ) )  ){
          index.t = (feature[[cur.ch]][,1] <= cur.end & feature[[cur.ch]][,1]>=cur.beg)
          temp.mean=c( temp.mean,  mean( feature[[cur.ch]][index.t, 2], na.rm=T)  ) 
          
        } 
        
        
      } #end loop through index
      
      if (!is.null(temp.mean) ){
        if (length(temp.mean)>1 ){
          cv_min = c( cv_min, sd(temp.mean, na.rm=T )/mean(temp.mean, na.rm=T) )
        } else{
          cv_min = c( cv_min, NA )
        }
        
        
      } else {
        cv_min = c( cv_min, NA )
      }
      
      
      
      
      cur.beg = cur.end #shift beginning of time window to end of last window
      
      
    } # end of loop through minutes
    
    
    #cv_time is average cv well
    cv_network =  c( cv_network, mean( cv_min, na.rm=T ) ) 
    
    
  }#end of loop through wells
  
  
  results = list( cv_time, cv_network)
  results
  
}  # end of .cv.feature function


#this function takes a spike object that has s[[1]]$allb, s[[1]]$bs and makes vector
# of events associated with each burst
# input: 1. s   2. feature = "IBI" or 'durn' or 'mean.isis'
# output: an object burst.train[[1]]:burst.train[[nChannels]] with timestamps at
# the beginning of each burst
.make.feature.vector<-function(s, feature, t.p ){
  
  if ( feature=="IBI" ){
    col.num=3
  } else if( feature =="durn" ){
    col.num=5
  }else if (feature == "mean.isis"){
    col.num = 6
  } else {
    print("No Feature selected: problem!")
  }
  
  
  #make a burst train
  feature.vector = list()
  feature.vector[[ length(s[[t.p]]$cw) ]]=0
  
  
  for (cur.ch in c(1:length( s[[t.p]]$cw) ) ){
    
    
    time= c(); feature = c() 
    if (length(s[[t.p]]$allb[[cur.ch]][,1] )>0){
      
      for (burst.n in c(1: length(s[[t.p]]$allb[[cur.ch]][,1] ) ) ){
        #this if statement is because the first row of s[[1]]$allb[[cur.ch]][1,3]==NA
        if ( ! (col.num==3 & burst.n==1 ) ){
          

        firstSoB = s[[t.p]]$allb[[cur.ch]][burst.n,1]
        time = c(time , s[[t.p]]$spikes[[cur.ch]][firstSoB] )
        feature=c(feature, s[[t.p]]$allb[[cur.ch]][burst.n, col.num])
        }
        
      }
      temp_3 = cbind (time, feature)
      
      feature.vector[[cur.ch]]= temp_3
      
    } else {
      feature.vector[[cur.ch]] = NA
      
    }
    
  } #end of loop through all channels
  
  feature.vector
} #end of .make.feature.vector 




# purpose: get minute to minute firing rate, and sd. to be able to 
#  compare to DMSO
# INPUT: s 
# output: 
# 1. s[[t.p]]$min.freq.well := the Mean(firing rate/min) across all channels in well
# 2. s[[t.p]]$min.freq.well.sd := sd( Mean(firing rate/min) across all channels in well )
# 3. s[[t.p]]$min.freq.ch := vector for each ch, each element isfiring rate for that minute
# 4. s[[t.p]]$min.freq.ch.sd := sd( s[[t.p]]$min.freq.ch ), one per channel

.make.freq.min<-function( s ){
 
  
for (t.p in 1:length(s) ){
  # min.freq.ch is a vector for each channel of min to min firing rates
  min.freq.ch = c()   
  
  # min.freq.sd is sd of min.freq.ch
  min.freq.ch.sd = c() 

  min.freq.well = c() #firing rate for each well
  
  min.freq.well.sd = c() #firing rate of each well

  # loop through each well
  for (w in c(1:length(unique(s[[t.p]]$cw) ) ) ){
    
    #get indices of current well
    index=which(s[[t.p]]$cw==unique(s[[t.p]]$cw)[w])
    
    #convert seconds to minutes (m), and round up 
    total.m = round( (s[[t.p]]$rec.time[2] - s[[t.p]]$rec.time[1])/60, 0)  #total minutes   

    #loop through channels in current well

    temp.ch.mean = c()
    for ( cur.ch in index ){
      #loop through minutes in current channel
      temp.mean=c()
      cur.beg = s[[t.p]]$rec.time[1]
      for (cur.min in c(1:(total.m-1) ) ){
        cur.end = cur.beg + 60 #look at 60 second window of time
        
        #n.min is used to get mean
        n.min = 1
        if ( cur.min == total.m-1 ){
          n.min = (s[[t.p]]$rec.time[2]%%60)/60 #this is modulus operator so remainder
        }

        temp.mean=c( temp.mean, 
                     sum( s[[t.p]]$spikes[[cur.ch]] <= cur.end & s[[t.p]]$spikes[[cur.ch]]>=cur.beg)/n.min )
        cur.beg = cur.end #shift beginning of time window to end of last window
      } #end loop through minutes
      
      min.freq.ch[[cur.ch]] = temp.mean 
      
      # sd 
      min.freq.ch.sd[[cur.ch]] = sd(temp.mean, na.rm = T )
      
      temp.ch.mean = c( temp.ch.mean, mean(min.freq.ch[[cur.ch]], na.rm=T ) )
      
    } # end of loop through index, through current well
    
    
    temp.well.mean = mean( temp.ch.mean, na.rm = T )
    # min.freq.well average min.freq.ch per well
    min.freq.well  =  c(min.freq.well, temp.well.mean  )
    
    # min.freq.well sd mean( min.freq.ch ) per well
    min.freq.well.sd  =  c(min.freq.well.sd, sd(temp.ch.mean, na.rm=T  ) )
    
  } # end of loop through wells
  
  # names !
  names(min.freq.ch) = names( s[[t.p]]$meanfiringrate)   

  names(min.freq.ch.sd) = names( s[[t.p]]$meanfiringrate) 
  
  names( min.freq.well ) = unique(s[[t.p]]$cw)
  
  names( min.freq.well.sd ) = unique(s[[t.p]]$cw)

  s[[t.p]]$min.freq.ch = min.freq.ch
  s[[t.p]]$min.freq.ch.sd = min.freq.ch.sd
  s[[t.p]]$min.freq.well = min.freq.well
  s[[t.p]]$min.freq.well.sd = min.freq.well.sd
  
} # end of for( t.p in 1:length(s) )
# return s
s
  
}  # end of .make.freq.min
    
    

# digitize spikes
# Purpose: to take a spike train and make it into 0,1 with a given time bin
# inputs:
# 1. a s, spike object
# 2. time.bin, a length for the time bin in milliseconds
# 3. t.p = time.point, so s[[t.p]]
# 4. w = well, which well do you want operation performed on
#   may be a name e.g. "F8", or index of well ie. 48
# output
# dig.s, dig.s[[i]] is a list of length floor(total time (ms)/time.bin (ms))
#    with entries of 0 1 in bin according to whether at least one spike occurned

.digi.spikes<-function( s, time.bin, t.p, w ){

    #get indices of current well
    if (is.character(w)){
      index.w=which(s[[t.p]]$cw==w )  
      
    } else {
      index.w=which(s[[t.p]]$cw==unique(s[[t.p]]$cw)[w] )
      
    }
    # get which channels are wanted
    channels.wanted = s[[t.p]]$channels[index.w]  
  
    #convert seconds to ms
    start.ms=round(s[[t.p]]$rec.time[1]*10^3,digits=0)
    stop.ms=round(s[[t.p]]$rec.time[2]*10^3,digits=0)
    total.ms=stop.ms-start.ms
    
    # get total number of bins
    total.bins = ceiling( total.ms/time.bin )   
    
    #digitize spikes
    dig.s=list()
    dig.s[[length(index.w)]]=rep(0, total.bins )
    names(dig.s)=channels.wanted
       
    #digitize spikes: loop through channels
    for (cur.ch in channels.wanted){
      
      dig.s[[cur.ch]]=rep(0, total.bins) # create a digital spike train
      
      # get current spike train in ms
      cur.nspikes = s[[t.p]]$nspikes[[cur.ch ]]
      
      # get spikes with initial time removed in ms
      cur.spike.train = 10^3 *( s[[t.p]]$spikes[[ cur.ch ]] - s[[t.p]]$rec.time[1] )
      
      for (cur.spike.ind in 1:cur.nspikes ){
        # get bin for spike
        cur.bin = ceiling(cur.spike.train[cur.spike.ind]/time.bin )
        #get millisecond
        dig.s[[cur.ch]][ cur.bin ] = 1
      }
    }
    
  dig.s
    
} #end .digi.spikes


# purpose: given 2 spike trains, make a new spike train of same time.bin length and number of bins
# that contains 1 if >1 spike occured in both bins, 0 else
# input:
# 1. digi.s1, a vector (element of list) of digitized spike train, output of .digi.spikes
# 2. digi.s2, same as above
.make.joint.digi.s<-function(dig.s1, dig.s2 ){
  
  # under null hypothesis of independent of channels you'd expect to see a certain group spike rate
  y_temp = dig.s1 + dig.s2 # element wise sum
  
  y_dig = floor( y_temp/2 ) #1 given time.bin if spike occured in both channels, 0 else
  
  
  y_dig
  
}


# purpose: this function makes a list for of rates for each channel in dig.s
#inputs
# dig.s, list of see .digi.spikes
.make.rates.for.digi.spikes<-function(dig.s ){
  
  rate = c() ;
  for (cur.ch in 1:length(dig.s) ){
    
    rate = c(rate, sum(dig.s[[cur.ch]])/length(dig.s[[cur.ch]]) )
    
  }
  
  names(rate) = names(dig.s)  
  rate
  
}

# wald stat for test of whether a sample proportion from bin (n, p unknown)
# equal to one given 
# 1. p_0, null prop (not sample prop)
# 2. p_sample, sample proportion
# 3. n, number of trials in binomial

.wald.stat.1sample.prop<-function(p_0, p_sample, n ){
  
  denom = sqrt( (p_sample* (1 - p_sample) )/ n )
  temp = ( p_sample - p_0 )/ denom
  
  temp2 = dnorm(temp, 0,1) # get p-value but comparing to N(0,1)
  
  z.w = c(temp, temp2)
  names(z.w) = c("Wald test value", "p-value")
  z.w
}









# wald stat for test of whether a sample proportion from bin (n, p unknown)
# equal to one given 
# 1. nsucess1, for sample 1
# 2. nsucess2, for sample 2
# 3. n1, n2 number of trails for sample 1,2
# output: z.w, test value and p-value
.wald.stat.2sample.prop<-function(nsucess1, nsucess2, n1, n2 ){
  
  p1 = nsucess1/n1 ;  p2 = nsucess2/n2; p.all = (nsucess1+nsucess2)/(n1+n2)
  
  denom = sqrt( (p.all* (1 - p.all) )*( (1/n1)+(1/n2) ) )
  temp = ( p1 - p2 )/ denom
  
  temp2 = dnorm(temp, 0,1) # get p-value but comparing to N(0,1)
  
  z.w = c(temp, temp2)
  names(z.w) = c("Wald 2-sample test value", "p-value")
  z.w
}



# purpose: to take all intervals between spikes in a given spike train
#  and return a new spike train of same ISI distribution but different actual times
# inputs:
# 1. spike train, a vector of spike times
# outputs:
# 1. a new spike train, with same dist'n of ISI but different timestamps
# example: new.spikes<- .shuffle.isis(s[[1]]$spikes[[1]])

.shuffle.isis<-function(spikes){
  
  isi.b = isi(spikes) # get the isi of a spike train
  # take a random sample, of nspikes size, without replacement
  isi.b.rs = sample(isi.b, size = length(isi.b), replace = FALSE, prob = NULL)
  
  new.spikes = c(); cur.spike = spikes[1] #initialize current spike
  for (spike in c(1:length(spikes) ) ){
    #fill in spike time 
    new.spikes = c(new.spikes, cur.spike )
    # update to new spike time
    cur.spike = cur.spike + isi.b.rs[spike] 
  }
  #return new spikes
  new.spikes
}

# purpose: compute network spikes for all channels 
# (active electrodes and non-active electrodes) alike
.plate.ns.ont<-function(s, t.p, ns.T, ns.N, sur ){
  
  ns=list()
  ns[[ length(unique( s[[t.p]]$cw ) ) ]] = c()
  for ( cur.well in unique( s[[t.p]]$cw ) ){
    
    ns[[cur.well]] <- .compute.ns(s[[t.p]], ns.T=ns.T, ns.N=ns.N, 
                                 sur=sur, whichcells=cur.well  )
    ns[[cur.well]]$well = cur.well
    print (paste("well ", cur.well, " finished" )  )
    
  } #end of for loop through wells in a plate
  ns
  
}







##
# purpose: to compute local correlation in ontogeny data
.local.corr.all.ont<-function(s, t.p) {
  
  # store data
  
  well.list<-unique(s[[t.p]]$cw)
  
  temp=list()
  #r.vector is a list for each channel
  r.vector=list()
  r.vector[[length( well.list )]]=0
  
  names(r.vector)= well.list
  
  #loop through wells
  for (cur.well in well.list ){
    
    #get indices of channels within current well
    index.ch=which(s[[t.p]]$cw==cur.well)
    
    #convert seconds to ms
    start.ms=round(s[[t.p]]$rec.time[1]*10^2,digits=0)
    stop.ms=round(s[[t.p]]$rec.time[2]*10^2,digits=0)
    total.ms=stop.ms-start.ms
    
    #digitize spikes
    dig.s=list()
    dig.s[[length(index.ch)]]=rep(0, total.ms )      
    
    #digitize spikes: loop through current channels
    #dig.s[[i]] is a list 1:# channels within current well
    for (cur.index in c(1:length(index.ch)) ){
      
      cur.ch = index.ch[cur.index]
      
      dig.s[[cur.index]]=rep(0, total.ms)
      
      current.s.ms=round(s[[t.p]]$spikes[[cur.ch]]*10^2, digits=0)- start.ms 
      for (j in current.s.ms ){
        #get millisecond
        dig.s[[cur.index]][ j ] = 1
      }
    } #end loop through current channels 
    
    
    #initialize c.vector: vector of mean(r), for each channel in current well
    
    
    #loop through channels in current well, getting cor between an electrode
    #and electrodes above, below, to left and to right
    # if those electrodes have AE
    c.vector = c(); ch.names=c();
    for (cur.ch in c(1:length(index.ch) ) ){
      
      temp.corr= c()
      for (other.ch in setdiff(c(1:length(index.ch )),cur.ch)  ){
        temp.corr=c(temp.corr, cor(dig.s[[cur.ch]], dig.s[[other.ch]]) )
        
      } #end of loop through other channels       
      
      #vector of mean correlations inside well by channel
      c.vector = c(c.vector, mean(temp.corr, na.rm=T ) )
      
    }#end of loop through channels
    
    names(c.vector)= s[[t.p]]$channels[index.ch]
    # each well get one vector
    r.vector[[cur.well]] = c.vector
    
  }#end of loop through wells
  
  r.vector   
  
}

# transfer entropy based off Liam Paninsky's suggestion on a paper by Ito
# it's basically the dynamic mutual entropy-is
# lag is the time lag at which you judge one against another
# code from users.utu.fi/attenka/TEpresentation081128.pdf
# .transfer.ent(X,Y)=0.011 means system X adds 0.011 digits of predictability to Y
#  .transfer.ent(Y,X)=0.044 means Y adds 0.044 digits of predictability to system X
# Y=c(0,0,1,1,0,1,1,1,0,0), X=(c(0,0,0,1,0,1,0,1,1,0)) then .transfer.ent(Y,X)
.transfer.ent<-function(Y,X,lag=1){
  
  L4=L1=length(X)-lag
  L3=L2=length(X)
  ##########################
  #1 p(Xn+s, Xn, Yn)
  ##########################
  TPvector1=rep(0,L1)
  for(i in 1:L1){
    TPvector1[i] = paste(c(X[i+lag], "i", X[i], "i",Y[i]), collapse="" )
  }
  
  TPvector1T=table(TPvector1)/length(TPvector1)
  
  ###############
  #2. p(Xn)
  ###############
  TPvector2=X
  TPvector2T=table(X)/sum(table(X))
  ############
  # 3. p(Xn,Yn)
  ###############
  TPvector3=rep(0,L3)
  for(i in 1:L3){
    TPvector3[i]=paste(c(X[i],"i",Y[i]), collapse="")
    
  }
  TPvector3T=table(TPvector3)/length(TPvector2)
  
  # 4 p(Xn+s, Xn)
  TPvector4=rep(0,L4)
  for(i in 1:L4){
    TPvector4[i]=paste(c(X[i+lag],"i",X[i]),collapse="")
  }
  
  TPvector4T=table(TPvector4)/length(TPvector4)
  
  #++++++++++++++++++++++++++
  # transfer entropy T(Y->X)
  #+++++++++++++++++++++++++++++
  SUMvector=rep(0,length(TPvector1T) )
  
  for(n in 1:length(TPvector1T) )
  {
    SUMvector[n] = TPvector1T[n] * log10((TPvector1T[n] * TPvector2T[(unlist(strsplit(names(TPvector1T)[n],
                                                                                      "i")))[2] ])/(TPvector3T[paste((unlist(strsplit(names(TPvector1T)[n],"i")))[2],"i",(
                                                                                        unlist(strsplit(names(TPvector1T)[n],"i")))[3],  sep="",  collapse="") ] * TPvector4T
                                                                                        [paste((unlist(strsplit(names(TPvector1T)[n],"i")))[1], "i", (unlist(strsplit(names(
                                                                                          TPvector1T)[n], "i" )))[2], sep="", collapse="")]))
    
  }
  return( sum(SUMvector) )
  
  
} #end of the .transfer.ent-function

