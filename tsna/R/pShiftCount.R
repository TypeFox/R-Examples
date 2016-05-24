# compute Gibson p-shift counts using carter's relevent package


pShiftCount<-function(nd,start=NULL,end=NULL, output=c('final','full')){
  if(!is.networkDynamic(nd)){
    stop('pShiftCount is only appropriate for dynamic networks')
  }
  if(!is.directed(nd)){
    stop('pShiftCount is only appropriate for directed networks')
  }
  output<-match.arg(output)
  requireNamespace('relevent')
  tel<-as.data.frame(nd,start=start,end=end)
  # sort tel into appropriate temporal order
  tel<-tel[order(tel[,1],tel[,2],tel[,3],tel[,4]),]
  # calculate a 'group' target
  # defining as any simulataneous events from the same source
  # the accum.ps function uses NA to indicate group target
  dupes <-which(duplicated(tel[c('onset','terminus','tail')]))
  group<-tel[,'head']
  group[dupes]<-NA
  group[dupes-1]<-NA  # the first element of the group was not flagged by duplicated function
  # append the group flag back onto tel
  tel<-cbind(tel,group)
  # now it is like a list of events to pass to relevent code
  events<-as.matrix(tel[c('onset','tail','group')])
  if(nrow(events)>1){
    counts<-relevent::accum.ps(events)
  } else {
    # make a dummy array
    counts<-matrix(NA,ncol=13,nrow=0,dimnames=list(list(),list("AB-BA","AB-B0","AB-BY","A0-X0","A0-XA","A0-XY","AB-X0","AB-XA","AB-XB","AB-XY","A0-AY","AB-A0","AB-AY")))
  }
  if(output=='final'){
    return(counts[nrow(counts),,drop=FALSE])
  } else {
    # return full matrix with pshift count for each event
    # append on the information about the spell and edge involved
    # for some reason counts has an additional row blank row at the beginning which we trim off
    return(cbind(counts[-1,],tel[,c('onset','terminus','tail','head')],group=is.na(group)))
    #return(counts) 
  }
}