setGeneric('unUsedRecords<-', function(obj, value){standardGeneric('unUsedRecords<-')})
setMethod('unUsedRecords<-', c(obj='.MoveTrackSingle', value='logical'), function(obj, value){
	  if(n.locs(obj)!=length(value))
		  stop('Selection length does not match with number of locations')
	  unUsed<-as(obj, '.unUsedRecords')
	  xNew<-obj[!value,]
	  xOld<-obj[value,]
	  df1<-unUsed@dataUnUsedRecords 
	  df2<-xOld@data
	  if(nrow(df1)!=0){
	  df2[,setdiff(names(df1),names(df2))] <- NA
	  df1[,setdiff(names(df2),names(df1))] <- NA
	  df3 <- rbind(df1,df2) 
	  }else{
		  df3<-df2
	  }
	  unUsedNew<-new('.unUsedRecords', 
		      timestampsUnUsedRecords=ifelse(is.null(unUsed@timestampsUnUsedRecords), list(xOld@timestamps),list(c(unUsed@timestampsUnUsedRecords, xOld@timestamps)))[[1]],   
		      sensorUnUsedRecords=factor(c(as.character(unUsed@sensorUnUsedRecords), as.character(xOld@sensor)), levels=levels(obj@sensor)),
		      dataUnUsedRecords=df3
		      ) 
	  new(class(obj), unUsedNew, xNew)
})
setMethod('unUsedRecords<-', c(obj='.MoveTrackStack', value='logical'), function(obj, value){
	  if(sum(n.locs(obj))!=length(value))
		  stop('Selection length does not match with number of locations')
	  unUsed<-as(obj, '.unUsedRecordsStack')
	  xNew<-obj[!value,]
	  xOld<-obj[value,]
	  df1<-unUsed@dataUnUsedRecords 
	  df2<-xOld@data
	  if(nrow(df1)!=0){
	  df2[,setdiff(names(df1),names(df2))] <- NA
	  df1[,setdiff(names(df2),names(df1))] <- NA
	  df3 <- rbind(df1,df2) 
	  }else{
		  df3<-df2
	  }
	  ts<-ifelse(is.null(timestamps(unUsed)), list(timestamps(xOld)),list(c(timestamps(unUsed), timestamps(xOld))))[[1]]
	  id<-factor(c(as.character(trackId(unUsed)), as.character(trackId(xOld))))
	  o<-order(id,ts)
	  unUsedNew<-new('.unUsedRecordsStack', 
		      timestampsUnUsedRecords=ts[o],   
		      sensorUnUsedRecords=factor(c(as.character(unUsed@sensorUnUsedRecords), as.character(xOld@sensor)))[o],
		      trackIdUnUsedRecords=id[o],
		      dataUnUsedRecords=df3[o,]
		      ) 
	  new(class(obj),  unUsedNew, xNew)
})

setGeneric('unUsedRecords', function(obj,...){standardGeneric('unUsedRecords')})
setMethod('unUsedRecords', c(obj='.unUsedRecordsStack'), function(obj, ...){
	  as(obj, '.unUsedRecordsStack')
})
setMethod('unUsedRecords', c(obj='.unUsedRecords'), function(obj, ...){
	  as(obj, '.unUsedRecords')
})

