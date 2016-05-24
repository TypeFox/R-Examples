#' Create a SpatialPointsDataFrame from an ejObject
#' 
#' @param x An ejObject
#' @export
as.SpatialPointsDataFrame.ejObject <- function(x){
	#convert all the records to spatial points dataframes
	recordList <- lapply(x$records, function(record){
				as.SpatialPointsDataFrame.ejRecord(record)
			})
	#bind the spatial and non-spatial data separately
	recordLocations <- do.call(rbind, lapply(recordList, function(x){x@coords}))
	recordData <- do.call(plyr::rbind.fill, lapply(recordList, function(x){x@data}))
	
	#return the object
	sp::SpatialPointsDataFrame(recordLocations, recordData)
}

#' Convert an ej attribute to a dataframe
#' 
#' @param x An ejAttribute
as.data.frame.ejAttribute <- function(x){
	result <- data.frame(x$value)
	names(result) <- x$name
	return(result)
}

#'convert an ejEvent to a SpatialPointsDataFrame
#' 
#' @param x An ejEvent object
as.SpatialPointsDataFrame.ejEvent <- function(x){
	#grab the columns
	eventCols <- do.call(cbind, c(list(eventId=x$id, name=x$name, date=x$date),lapply(x$attributes, as.data.frame.ejAttribute)))
	sp::SpatialPointsDataFrame(x$location, eventCols)
}

#'convert an ejRecord to a SpatialPointsDataFrame
#' 
#' @param x an ejRecord object
as.SpatialPointsDataFrame.ejRecord <- function(x){
	#convert all the events to SpatialPointsDataFrames
	eventList <- lapply(x$events, function(event){
		#only work on those events with a location
		result <- NULL
		if(class(event$location) == "SpatialPoints"){
			result <- as.SpatialPointsDataFrame.ejEvent(event)
		}
		result
	})
	
	eventLocations <- do.call(rbind, lapply(eventList, function(x){x@coords}))
	eventData <- do.call(plyr::rbind.fill, lapply(eventList, function(x){x@data}))

	#convert the indicidual attributes to columns
	recordAttributes <- do.call(cbind, c(list(recordId=x$id),lapply(x$attributes, as.data.frame.ejAttribute)))
	
	
	#here we do something a bit different. Normally we one record per row
	#but for a spatial object we are going to repeat records for each
	#spatial event
	
	#duplicate to match the event data
	recordAttributes <- recordAttributes[rep(1,nrow(eventData)),]
	
	resultDF <- cbind(recordAttributes, eventData)
	#construct the Spatial points dataframe
	sp::SpatialPointsDataFrame(eventLocations, resultDF)
}