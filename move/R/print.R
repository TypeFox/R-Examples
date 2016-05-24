###Print function for a Move and MoveStack object
setGeneric("print")
setMethod("print", ".MoveTrack", function(x){
	  getMethod("print","Spatial")(x) 
	  timeRange <- range(x@timestamps)
	  cat("timestamps  :",paste(timeRange, collapse=" ... "),capture.output(round(difftime(timeRange[2],timeRange[1]))), " (start ... end, duration) \n")  
	  cat("sensors     :",paste(as.character(unique(x@sensor)), collapse=", "),"\n")
	  cat("indiv. data :", paste(collapse=", ", names(x@idData)),"\n")
})
setMethod("print", "dBMvarianceTmp", function(x){
	  cat("margin      :", x@margin,'\n')
	  cat("window size :", x@window.size,'\n')# maybe add means usw here later
})
setMethod("print", ".MoveTrackSingle", function(x){
	  callNextMethod()
	  cat("indiv. value:", paste(apply(x@idData, 2,min)),"\n")
	  print(unUsedRecords(x))
})
setMethod("print", ".MoveTrackStack", function(x){
	  callNextMethod()
	  cat("min ID Data :", paste(apply(x@idData, 2,min)),"\n")
	  cat("max ID Data :", paste(apply(x@idData, 2,max)),"\n")
	  cat("individuals :", paste(as.character((unique(x@trackId)[1:min(c(15,n.indiv(x)))])), collapse=", "),"\n")
	  print(unUsedRecords(x))
})
setMethod('print','.unUsedRecords', function(x){
	  try(silent=TRUE, if(length(x@timestampsUnUsedRecords)>=1)
	  cat("unused rec. :", length(x@timestampsUnUsedRecords),"\n") )
})
setMethod("print", "dBMvariance", function(x){
	  callNextMethod()
	  print(as(x,"dBMvarianceTmp"))
})
setMethod("print", "Move", function(x){
	  callNextMethod()
	  print(as(x,".MoveGeneral"))
})
setMethod("print", "MoveBurst", function(x){
	  callNextMethod()
	  print(as(x,".MoveGeneral"))
})
setMethod("print", "MoveStack", function(x){
	  callNextMethod()
	  print(as(x,".MoveGeneral"))
})
setMethod("print",".MoveTrackSingleBurst", function(x){
	  callNextMethod()
	  cat("bursts      :", paste(levels(x@burstId),': ',unlist(lapply(lapply(levels(x@burstId), '==', x@burstId), sum)), sep='', collapse=', '),"\n")
})
setMethod("print",".MoveGeneral", function(x){
	  if(length(x@license)!=0)
		  if(nchar(x@license)>(getOption('width')-18)){
		  cat("license     :",substr(x@license,1, getOption('width')-18),"...\n")
		  }else{
		  cat("license     :",x@license,"\n")
	  }
	  if(length(x@citation)!=0)
		  cat("citation     :",x@citation,"\n")
	  if(length(x@study)!=0)
		  cat("study name   :",x@study,"\n")
	  cat("date created:",format(x@dateCreation),"\n")

})
