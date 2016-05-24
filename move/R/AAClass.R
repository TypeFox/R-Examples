setClassUnion(".OptionalPOSIXct", c("POSIXct","NULL"))
lapply(lapply(c('url','file','gzfile','unz','pipe','fifo'),c,'connection'), setOldClass,   where = environment())
setOldClass("ltraj")
setClass(Class = ".MoveGeneral",
	 representation = representation(
					 dateCreation = "POSIXct",
					 study = "character",
					 citation = "character",
					 license = "character"),
	 prototype = prototype(
			       dateCreation = Sys.time(),
			       study = as.character(),
			       citation = as.character(),
			       license = as.character()),
	 validity = function(object){
		 if(length(object@citation)>1)
			 stop("Citation has length unequal to 0 or 1")
		 if(length(object@license)>1)
			 stop("License has length unequal to 0 or 1")
		 return(TRUE)
	 }
	 )
setClass(Class='.unUsedRecords',
	 representation=representation(
				       timestampsUnUsedRecords = ".OptionalPOSIXct",
				       sensorUnUsedRecords="factor",
				       dataUnUsedRecords='data.frame'
				       ),
	 prototype = prototype(
			       timestampsUnUsedRecords = NULL,
			       sensorUnUsedRecords=factor(),
			       dataUnUsedRecords=rbind(data.frame())# rbind needed so that the unit check thinks it is the same after deparsing
			       ),
	 validity = function(object){
		 if(length(object@sensorUnUsedRecords)!=nrow(object@dataUnUsedRecords))
			 stop("Sensors and data length for unused records do not match")
		 if(length(object@timestampsUnUsedRecords)!=nrow(object@dataUnUsedRecords))
			 stop("Timestamps and data length for unused records do not match")
		 return(TRUE)
	 }
	 )
setClass(Class='.unUsedRecordsStack',contains='.unUsedRecords',# maybe at somestage introduce checking of duplicate timestamps in unused records
	 representation=representation(
				       trackIdUnUsedRecords='factor'
				       ),
	 prototype = prototype(
			       trackIdUnUsedRecords=factor()
			       ),
	 validity = function(object){

		 if(length(object@trackIdUnUsedRecords)!=nrow(object@dataUnUsedRecords))
			 stop("TrackId and data length for unused records do not match")
		 if(length(object@timestampsUnUsedRecords)>0)
		 if(!all((diffs<-unlist(tapply(object@timestampsUnUsedRecords, list(object@trackIdUnUsedRecords, droplevels(object@sensorUnUsedRecords)),diff)))>0))
		 {
			 if(any(diffs==0)){
		   tmp<-duplicated(cbind(object@timestampsUnUsedRecords, object@trackIdUnUsedRecords, object@sensorUnUsedRecords))# made new test in check for higher speed but have this one still here for more clear reporting
			 stop("The data set includes double timestamps per ID in the unused records (first one:", object@trackIdUnUsedRecords[tmp][1]," ",object@sensorUnUsedRecords[tmp][1]," ",object@timestampsUnUsedRecords[tmp][1], ")")
		   }else{
			 stop('The data set includes un ordered timestamps in the unUsedRecordsStack')
			   }
		 }
		 return(TRUE)
	 }
	 )


setClass(Class = ".MoveTrack",contains=c("SpatialPointsDataFrame"),
	 representation = representation(
					 timestamps = "POSIXct",
					 idData = "data.frame",
					 sensor="factor"
					 ),
	 prototype = prototype(
			       timestamps = as.POSIXct(NA),
			       idData = data.frame(),
			       sensor=factor()
			       ),
	 validity = function(object){
		 if(length(object@timestamps)!=nrow(object@coords))
			 stop("Number of timestamps does not match the number of coordinates")
		 if(any(is.na(object@timestamps)))
			 stop("There are NA timestamps records")
		 if(length(object@sensor)!=nrow(object@coords))
			 stop("Number of sensors observations does not match the number of coordinates")
		 if(any(is.na(object@sensor)))
			 stop("There are NA sensor records")
		 return(TRUE)
	 }
	 )



setClass(Class = ".MoveTrackSingle",contains=c(".MoveTrack",'.unUsedRecords'), 
	 representation = representation(),
	 prototype = prototype(),
	 validity = function(object){
		 if(any((tmp <- timeLag(object))<0))
			 stop("The dataset includes unsorted timestamps")
		 tmp <- c(T, tmp==0)|c(tmp==0,T)#only check for those that have the same time if the sensor is different
		 if (any(dups <- duplicated(data.frame(format(object@timestamps[tmp],"%Y %m %d %H %M %OS4"), object@sensor[tmp]))))
			 stop("The dataset includes double timestamps first one:", object@timestamps[tmp][dups][1], ")")
		 if(nrow(idData(object, drop=F))>1)
			 stop("More than 1 row are stored in the idData data.frame")
		 if(nrow(idData(object, drop=F))<1)
			 stop("Less than 1 row are stored in the idData data.frame")
		 if(!identical(levels(object@sensorUnUsedRecords),levels(object@sensor)))
			 stop('Levels of unused records dont match with sensor')
		 timestampsUnUsedDuplicated<-object@timestampsUnUsedRecords[object@timestampsUnUsedRecords %in% object@timestamps]
		 if(length(timestampsUnUsedDuplicated)!=0)
		 {
		 s<-c(object@timestamps, object@timestampsUnUsedRecords)%in% timestampsUnUsedDuplicated
		 	if (any(dups <- duplicated(data.frame(format(c(object@timestamps, object@timestampsUnUsedRecords)[s],"%Y %m %d %H %M %OS4"), 
							       c(as.character(object@sensor), as.character(object@sensorUnUsedRecords))[s]))))
				 stop("A timestamp of an unused record coincides with a normal timestamp")

		 }
		 #this check cant work since coordinates columns are not present in data maybe look for solution
		# if(any(names(object@data)!=names(object@dataUnUsedRecords)))
		#	 stop('names of data and unused data records dont match')
		# if(any(unlist(lapply(object@data, class))!=unlist(lapply(object@dataUnUsedRecords, class))))
		#	 stop('classes for data and unusedrecords data dont match')
		 return(TRUE)
	 }
	 )


setClass(Class = "Move", contains=c(".MoveTrackSingle",".MoveGeneral"),
	 representation = representation (
					  ),
	 prototype = prototype(
			       ),
	 validity = function(object){
		 return(TRUE)
	 }
	 )


setClass(Class = ".MoveTrackStack", contains = c(".MoveTrack", ".unUsedRecordsStack"),
	 representation = representation(
					 trackId = "factor"),
	 prototype = prototype(
			       trackId = factor()),
	 validity = function(object){
		 if(length(object@trackId)!=nrow(object@coords))
			 stop("Length of trackId does not match the number of coordinates")
		 if(!all(unlist(tapply(object@timestamps, list(trackId(object), droplevels(object@sensor)),diff))>0))
		 {
			 tmp<-duplicated(cbind(object@timestamps, trackId(object), object@sensor))# made new test in check for higher speed but have this one still here for more clear reporting
			 stop("The data set includes double timestamps per ID (first one:", object@trackId[tmp][1]," ",object@sensor[tmp][1]," ",object@timestamps[tmp][1], ")")
		 }
		 if(any(unlist(lapply(tapply(object@timestamps,object@trackId, order),diff))!=1))
			 stop("Not ordered timestamps per individual occured")
		 if(any(levels(object@trackId)!=validNames(levels(object@trackId))))
			 stop('No good names for trackId levels')
		 if(length(unique(object@trackId))!=nrow(object@idData))
			 stop("Not same number of unique IDs and rows in the idData data.frame")
		 if(any(is.na(object@trackId)))
			 stop("There are NA trackId records")
		 if(any(sort(as.character(unique(object@trackId)))!=sort(unique(rownames(object@idData))))){
			 stop("No match between rownames in idData and ids along track")} 
		 if(!all(unique(object@trackIdUnUsedRecords)%in%unique(object@trackId)))
			 stop("There are records for individuals where no real records are present")
		 if(!identical(levels(object@sensorUnUsedRecords),levels(object@sensor)))
			 stop('Levels of unused records dont match with sensor')
		 if(!identical(levels(object@trackIdUnUsedRecords),levels(object@trackId)))
			 stop('Levels of unused records dont match with trackId')
		 timestampsUnUsedDuplicated<-object@timestampsUnUsedRecords[object@timestampsUnUsedRecords %in% object@timestamps]
		 if(length(timestampsUnUsedDuplicated)!=0)
		 {
		 s<-c(object@timestamps, object@timestampsUnUsedRecords)%in% timestampsUnUsedDuplicated
		 	if (any(dups <- duplicated(t<-cbind(
							  format(c(object@timestamps, object@timestampsUnUsedRecords)[s],"%Y %m %d %H %M %OS4"), 
							  c(as.character(object@sensor), as.character(object@sensorUnUsedRecords))[s], 
							  c(as.character(trackId(object)),as.character( object@trackIdUnUsedRecords))[s]))))
				 stop("A timestamps of a unused record coincides with a normal timestamps")

		 }
		 if(sum(diff(as.numeric(trackId(object)))!=0)!=(nrow(idData(object, drop=F))-1))
			 stop('The data in the MoveStack object are not grouped per individual')
		 if(any(as.character(unique(trackId(object)))!= rownames(idData(object, drop=F))))
			 stop('Order of objects in the idData is not the same as in the trackId')
		 if(!identical(as.character(unique(trackId(object))), levels(object@trackId)))
			 stop('Order of levels in the trackId should be same as order of individuals') 
		 #this check cant work since coordinates columns are not present in data
		# if(any(names(object@data)!=names(object@dataUnUsedRecords)))
		#	 stop('names of data and unused data records dont match')
		# if(any(unlist(lapply(object@data, class))!=unlist(lapply(object@dataUnUsedRecords, class))))
		#	 stop('classes for data and unusedrecords data dont match')
		 return(TRUE)
	 }
	 )



setClass(Class = "MoveStack", contains = c(".MoveTrackStack",".MoveGeneral"),
	 representation = representation(),
	 prototype = prototype(),
	 validity = function(object){
		 return(TRUE)
	 }
	 )
setClass(Class = ".MoveTrackSingleBurst", contains = c(".MoveTrackSingle"), 
	 representation = representation(
					 burstId = "factor"), 
	 prototype = prototype(
			       burstId = factor()), 
	 validity = function(object) {
		 if(any(levels(object@burstId)!=validNames(levels(object@burstId))))
			 stop('no good names')
		 if(length(object@burstId)!=(length(object@timestamps)-1))
			 stop("Burst ids need to be one shorter than the number of coordinates since it is a segment property")
		 return(TRUE)
	 })


setClass(Class = "MoveBurst", contains = c(".MoveTrackSingleBurst", ".MoveGeneral"), 
	 validity = function(object) {
		 return(TRUE)
	 })


setClass(Class = "dBMvarianceTmp", 
	 representation = representation(
					 window.size = "numeric", 
					 margin = "numeric", 
					 means = "numeric", 
					 in.windows = "numeric", 
					 interest = "logical", 
					 break.list = "numeric"), 
	 prototype = prototype(
			       window.size = numeric(), 
			       margin = numeric(), 
			       means = numeric(), 
			       in.windows = numeric(), 
			       interest = logical(), 
			       break.list = numeric()), 
	 validity = function(object) {
		 if (length(unique(lengths<-c(length(object@means), length(object@in.windows), length(object@interest)))) != 1) 
			 stop("Length does not match between means, in.windows and interest (", paste(lengths, collapse=', '),')')
		 if (length(object@margin) != 1) 
			 stop("Margin length not 1")
		 if (length(object@window.size) != 1) 
			 stop("Window size length not 1")
		 if( any(is.na(object@means[object@interest])))
			 stop('There are not variance estimates for segments of interest')
		 if( tail(object@interest,1))
			 stop('The last value of interest cant be true since interest has the same length as the number of locations in the object and the last value thus does not refer to a segment')
		 return(TRUE)
	 })


setClass(Class = "dBMvariance", contains = c(".MoveTrackSingle", "dBMvarianceTmp"), 
	 validity = function(object) {
		 if (length(object@means) != nrow(object@coords)) 
			 stop("Number of coordinates does not match the number of means")
		 return(TRUE)
	 })

setClass(Class = "dBMvarianceBurst", contains = c(".MoveTrackSingleBurst", "dBMvarianceTmp"), 
	 validity = function(object) {
		 if (length(object@means) != nrow(object@coords)) 
			 stop("Number of coordinates does not match the number of means")
		 return(TRUE)
	 })

setClass(Class = "dBMvarianceStack", contains = c(".MoveTrackStack", "dBMvarianceTmp"), 
	 validity = function(object) {
		 if (length(object@means) != nrow(object@coords)) 
			 stop("Number of coordinates does not match the number of means")
		 return(TRUE)
	 })
########################
########## UD ##########
########################


setClass(Class = ".UDStack", contains = c("RasterStack"), 
	 representation = representation(method = "character"), 
	 prototype = prototype(
			       method = as.character()), 
	 validity = function(object) {
		 #if (!all(apply(values(object), MARGIN = 2, FUN = function(X) isTRUE(all.equal(sum(X), 1, check.attributes=F))))) 
		 if(!all.equal(rep(1,nlayers(object)),cellStats(object, sum), check.attributes=F)) 
			 stop("One or more of the used rasters are not a UD, because they sum not to 1)")
	 })

setClass(Class = ".UDBurstStack", contains = c("RasterStack"), 
	 representation = representation(method = "character"), 
	 prototype = prototype(
			       method = as.character()), 
	 validity = function(object) {
		 #if (!all(apply(values(object), MARGIN = 2, FUN = function(X) isTRUE(all.equal(sum(X), 1, check.attributes=F))))) 
		 if(!all.equal(1,sum(s<-cellStats(object,'sum')), check.attributes=F)) 
			 stop("All values in a burst stack need to sum up to one")
		 if(class(z<-getZ(object))!='difftime')
			 stop('The Z vector needs to represent the time contribution in the for of a difftime vector')
		 if(!
		    all.equal( s , as.numeric(z, units='mins')/sum(as.numeric(z, units='mins')), check.attributes=F)
		    )
			 stop('The Z vector need to correspond to the sum of the layers')
		 return(TRUE)
	 })

setClass(Class = ".UD", contains = c("RasterLayer"), 
	 representation = representation(method = "character"), 
	 prototype = prototype(
			       method = as.character()), 
	 validity = function(object) {
		 if (!isTRUE(all.equal(tmp<-sum(values((object))), 1))) 
			 stop("The used raster is not a UD (sum unequal to 1), sum is: ", sprintf("%.15f",tmp)," One possible cause is loss of accuracy due to writing raster to disk with dataType FLT4S this can be solved preventing disk usage or changing data type")
		 return(TRUE)
	 })


### Defining the class of the Brownian Bridge Movement Model object
setClass(Class = "DBBMMStack", contains = c(".UDStack"), 
	 representation = representation(
					 DBMvar = "dBMvarianceStack", 
					 ext = "numeric"), 
	 prototype = prototype(
			       ext = as.numeric()), 
	 validity = function(object) {
		 if (!all(unique(object@DBMvar@trackId) == names(object))) 
			 stop("The layer names of the raster objects do not match the trackIDs of the DBMvarStack")
	 })


setClass(Class = "DBBMM", contains = c(".UD"), 
	 representation = representation(
					 DBMvar = "dBMvariance", 
					 ext = "numeric"), 
	 prototype = prototype(
			       ext = as.numeric())
	 )
setClass(Class = "DBBMMBurstStack", contains = ".UDBurstStack", 
	 representation = representation(
						 DBMvar = "dBMvarianceBurst", 
						 ext = "numeric"), 
	 prototype = prototype(ext = as.numeric()), 
	 validity = function(object) {
		 validObject(object@DBMvar)
		 return(T)
	 })
## dynBGB classes ##
setClass(Class = "dBGBvarianceTmp", 
	 representation = representation(windowSize = "numeric", 
					 margin = "numeric", paraSd = "numeric", orthSd = "numeric", nEstim = "numeric", 
					 segInterest = "logical"), prototype = prototype(windowSize = numeric(), margin = numeric(), 
					 paraSd = numeric(), orthSd = numeric(), nEstim = numeric(), segInterest = logical()), 
	 validity = function(object) {
		 if (length(unique(c(length(object@paraSd), length(object@orthSd), length(object@nEstim), 
				     length(object@segInterest)))) != 1) 
			 stop("Length is not the same between segInterest, orthSd, paraSt, nEstim")
		 if (length(object@margin) != 1) 
			 stop("Margin length not 1")
		 if (length(object@windowSize) != 1) 
			 stop("Window size length not 1")
		 if (2 != sum(uneq <- rev(object@segInterest) != object@segInterest)) 
			 stop("something wrong with segInterest")
		 if (any(uneq != rev(uneq))) 
			 stop("something wrong with segInterest")
		 return(TRUE)
	 })

setClass(Class = "dBGBvariance", contains = c(".MoveTrackSingle", "dBGBvarianceTmp"), 
	 validity = function(object) {
		 if (length(object@segInterest) != nrow(object@coords)) 
			 stop("Number of coordinates does not match the number of means")
		 return(TRUE)
	 })
setClass("dynBGB", contains=c(".UD"), representation(var="dBGBvariance"))

