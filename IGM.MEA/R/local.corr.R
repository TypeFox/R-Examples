#this function takes in a list of spikes (s), with burst info
# digitizes spikes into 10ms bins (1 in bin if any spike, 0 otherwise)
# then computes correlation btw ae above, below, to left and to right
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
