
#' This function takes a single record from the obkClass data and 
#' converts to the epiJSON format
#' 
#' @param x An record from the obkData 
#' @examples
#' \dontrun{
#' #because this function is not exported this example won't work outside the package
#' require('OutbreakTools')
#' data(ToyOutbreak)
#' x <- subset(ToyOutbreak,1)  
#' processRecord(x)
#' }
#' 
#' @return an ejRecord
processRecord <- function(x){
	#get the record ID
	recordID <- row.names(x@individuals)
		
	#convert to attributes
	attributes <- dataFrameToAttributes(x@individuals)
	
	#process the record frames to events
	recordFrames <- x@records
	events<-c()
	eventID <- 1
	for(recordFrame in names(recordFrames)){
		#skip empty frames
		if (nrow(recordFrames[[recordFrame]])!=0){
			events <- c(events, processRecordFrame(recordFrames[[recordFrame]], recordFrame, eventID))
			eventID <- eventID + 1
		}
	}
	#fix the event ids
	#events <- lapply(seq_along(events), function(i){x<-events[[i]]; x$id <- i; x})
	create_ejRecord(id=recordID, attributes, events)
}


#' This function processes events from the obkClass data and 
#' converts to the epiJSON format
#' 
#' @param x An record from the obkData 
#' @param recordFrameName The event of interest
#' @param eventID The id of the event
#' @examples
#' \dontrun{
#' #because this function is not exported this example won't work outside the package
#' require('OutbreakTools')
#'
#' data(ToyOutbreak)
#' x=subset(ToyOutbreak,2)@@records[[1]]
#' processeventFrame(x,"Fever")
#' }
#' @return an ejEvent
processRecordFrame <- function(x, recordFrameName, eventID){	
	lapply(1:nrow(x), function(i){
		eventAttributes <- dataFrameToAttributes(x[i,3:ncol(x), drop=FALSE])
		create_ejEvent(id=eventID, date=as.POSIXct(x$date[i]), name=recordFrameName, attributes=eventAttributes)
	})	
}


#' This function processes objects from the obkClass data and 
#' converts to the epiJSON format
#' 
#' @param x An record from the obkData 
#' @param metadata The list of the components in the metadata
#' @param ... other parameters (to maintain consistency with the generic)
#' @note There is a slight mismatch in symantics here obkData individuals are
#'  equivelent to EpiJSON records and obkData records are EpiJSON events. This
#'  is beause in EpiJSON the unit of record is not necessarily an individual
#'  (it could, for example, be a region or hospital, etc). 
#' @examples
#' require('OutbreakTools')
#' data(ToyOutbreak)
#' x=subset(ToyOutbreak,2)
#' as.ejObject(x, metadata=list())
#' 
#' @return an ejObject
#' @method as.ejObject obkData
#' @export 
as.ejObject.obkData <- function(x, metadata=list(), ...){
	records <- lapply(OutbreakTools::get.individuals(x), function(xx){
		processRecord(OutbreakTools::subset(x, xx))		
	})

	create_ejObject(metadata=metadata, records=records)
}

# #' Create an obkData object from an ejObject
# #' 
# #' @export 
# as.obkData.ejObject <-function(){}