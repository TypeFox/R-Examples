#require(move)
setGeneric("interpolateTime", function(x, time,spaceMethod=c('euclidean','greatcircle','rhumbline'),...){standardGeneric("interpolateTime")})
setOldClass("difftime") 
setMethod("interpolateTime", 
          signature=c(".MoveTrackSingle",'numeric'), 
          definition=function(x,time, ...){
		  stopifnot(length(time)==1)
		time<-seq(min(timestamps(x)), max(timestamps(x)),length.out=time)
		callGeneric()
	  })
setMethod("interpolateTime", 
          signature=c(".MoveTrackSingle",'difftime'), 
          definition=function(x,time, ...){
		  stopifnot(length(time)==1)
		time<-seq(min(timestamps(x)), max(timestamps(x)),by=time)
		callGeneric()
	  })
setMethod("interpolateTime", 
          signature=c(".MoveTrackSingle",'POSIXct'), 
          definition=function(x,time, spaceMethod=c('euclidean','greatcircle','rhumbline'),...){
		  stopifnot(max(time)<=max(timestamps(x)))
		  stopifnot(min(time)>=min(timestamps(x)))
		  spaceMethod<-match.arg(spaceMethod)
		  if(spaceMethod=='euclidean' & isLonLat(x))
			  warning('Euclidean interpolation seems unsuitable for the longitude latitude projection')
		  # the next two steps are the slow ones
		  #previous loc
		  prevLoc<-unlist(lapply(lapply(lapply(time, '>=', timestamps(x)), which),max))
		  # next loc
		  nextLoc<-unlist(lapply(lapply(lapply(time, '<=', timestamps(x)), which),min))
		  p<-as.numeric(time-timestamps(x)[prevLoc], units='secs')/ as.numeric(timestamps(x)[nextLoc]-timestamps(x)[prevLoc],units='secs')
		  if(all(is.nan(p)))
			  return(x[prevLoc,])
		  fun<-switch(spaceMethod,
			      euclidean=function(x,y,p){(x)*(1-p)+(y)*p},
			      greatcircle=function(x,y,p){destPoint(x,ifelse(is.na(bearing(x, y)),0,bearing(x,y)),distHaversine(x, y)*p)},
			      rhumbline=function(x,y,p){ destPointRhumb(x,ifelse(is.na(bearingRhumb(x, y)),0,bearingRhumb(x,y)),distRhumb(x, y)*p) }
			      )
	#	  crds<-do.call('rbind',mapply(function(x,y,p,f,m){if(is.nan(p)){return(coordinates(m[x,]))}else{f(m[x,],m[y,],p)}},x=prevLoc, y=nextLoc,p=p,MoreArgs=list(m=x,f=fun), SIMPLIFY=FALSE))

		  crds<-matrix(fun(coordinates(x)[prevLoc,], coordinates(x)[nextLoc,], ifelse(is.nan(p),0, p)),ncol=2)# matrix incase only one location
	colnames(crds)<-colnames(coordinates(x))
  sensor<-factor(ifelse(is.nan(p),as.character(x@sensor[prevLoc]),'interpolateTime'), levels=c(levels(x@sensor),'interpolateTime'))
		  m<-new('Move', x[prevLoc,],coords=crds, timestamps=time, sensor=sensor, sensorUnUsedRecords=factor(x@sensorUnUsedRecords, levels=c(levels(x@sensorUnUsedRecords),'interpolateTime')) )
		  m[!is.nan(p),]<-NA
		  return(m)
          })

