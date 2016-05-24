#' Convert a dataframe to an ejObject
#'
#'
#' @param x The dataframe to convert
#' @param recordID an ID for the record, if NA one is created
#' @param recordAttributes A character vector containing the names of the
#'  columns in the dataframe that are attributes of the record
#' @param eventDefinitions A list of event definitions
#' @param metadata A list of metadata ejAttribute objects describing the dataset
#' @param ... other parameters (to maintain consistency with the generic)
#' @note We assume one row per record.
#' @examples 
#' #Here we use the toyll (toy line-list) data provided inside the package
#' data(toyll)
#' toyll
#' #Here we some data clean-up (make the date columns POSIX objects)
#' #Date-time conversion could be made automatic but introduces
#' #a big risk of mis-conversion so we leave this to the user
#' #to ensure that we don't silently corrupt their data
#' 
#' toyll$date.of.onset <- as.POSIXct(toyll$date.of.onset)
#' toyll$date.of.admission <- as.POSIXct(toyll$date.of.admission)
#' toyll$date.of.discharge <- as.POSIXct(toyll$date.of.discharge)
#' toyll$contact1.date <- as.POSIXct(toyll$contact1.date)
#' toyll$contact2.date <- as.POSIXct(toyll$contact2.date)
#' toyll$contact3.date <- as.POSIXct(toyll$contact3.date)
#' 
#' ind.fields <- c(names(toyll)[1:5], "hospital", "fever", "sleepy")
#' x <- as.ejObject(toyll,
#'                  recordAttributes=ind.fields,
#'                  eventDefinitions=list(
#'                  define_ejEvent(name="admission", date="date.of.admission"),
#'                  define_ejEvent(name="discharge", date="date.of.discharge"),
#'                  define_ejEvent(name="contact1", date="contact1.date", attributes="contact1.id"),
#'                  define_ejEvent(name="contact2", date="contact2.date", attributes="contact2.id"),
#'                  define_ejEvent(name="contact3", date="contact3.date", attributes="contact3.id")
#'                  ))
#' print(x)
#' @export
#'
as.ejObject.data.frame <- function(x, recordID=NA, recordAttributes, eventDefinitions, metadata=list(), ...){
	#iterate over the dataframe and create an record event for each row
	records <- lapply(1:nrow(x), function(i){
				#grab the attributes
				attributeDF <- x[i,unlist(recordAttributes), drop=FALSE]
				if (!is.null(names(recordAttributes))){
					newNames <- names(attributeDF)
					newNames[names(recordAttributes) != ""] <- names(recordAttributes)[names(recordAttributes) != ""]
					names(attributeDF) <- newNames
				}
				
				attributes <- dataFrameToAttributes(attributeDF)
				
				#now work over the eventDefinitions to get the events
				events <- lapply(1:length(eventDefinitions), function(j){
							rd <- eventDefinitions[[j]]
							
							#for the event attributes
							eventAttributeDF <- x[i,unlist(rd$attributes), drop=FALSE]
							if (!is.null(names(rd$attributes))){
								newNames <- names(eventAttributeDF)
								newNames[names(rd$attributes) != ""] <- names(rd$attributes)[names(rd$attributes) != ""]
								names(eventAttributeDF) <- newNames
							}
							eventAttributes <- dataFrameToAttributes(eventAttributeDF)
							
							#generate location
							location <- NULL
							if(!is.null(rd$location)){
								location <- sp::SpatialPoints(x[i, unlist(rd$location[c("x","y"), drop=FALSE])], proj4string=sp::CRS(rd$location$proj4string))
							}
							#generate te date
							date <- notNA(x[i, rd$date], x[i, rd$date], NULL)
							
							#if we have a date or location then create the event
							if (any(!is.null(c(date,location)))){
								create_ejEvent(
										id=notNA(rd$id, x[i,rd$id], j),
										name=notNA(x[i, rd$name],x[i, rd$name], rd$name),
										date=date,
										location=location,
										attributes=eventAttributes
								)
							} else {
								NULL
							}
						})
			
				#grab the record id
				id <- ifelse(is.na(recordID), i, x[i,recordID])
				
				#create and return the record
				create_ejRecord(id, attributes[!sapply(attributes, is.null)], events[!sapply(events, is.null)])
			})
	create_ejObject(metadata, records)
}

#' Creates a event definition
#'
#' Simplifies the definition of events from columns within a dataframe
#' @param id A character string naming the column that defines the id for a
#'  event. May be NA, and if so will be automatically generated.
#' @param name Either a character string with the event name. Or the name of the
#'  column that contains names for events as a character string. If the column 
#'  name is found in the input data.frame to \code{\link{as.ejObject.data.frame}} 
#'  then the data in the column is used otherwise the string itself is used.
#'  Event names might be things such as infection, swab, hospital admission,etc. 
#' @param date A character string naming the column that defines the date
#'  an event occured. This should be in POSIXct format. May be NA.
#' @param  location A list with entities x, y and proj4string. x and y should be
#'  character strings naming the columns where the x and y of the location are
#'  defined. crs may be "" or a proj4string.
#' @param attributes A character vector naming the columns for attributes of the
#'  event. The attributes will be named after the columns, with type taken from
#'  column type.
#' @export
define_ejEvent <- function(id=NA, name=NULL, date=NULL, location=NULL, attributes=NULL){
	if(anyNA(name) || is.null(name))
		stop("The name for an event must be supplied in define_ejEvent")
	if(is.null(location) && is.null(date))
		stop("An event requires either a date or a location")
	if(anyNA(location) && anyNA(date))
		stop("An event requires either a date or a location")
	structure(list(
					id=id,
					name=name,
					date=date,
					location=location,
					attributes=attributes
			), class="ejDFeventDef")
}


#' convert an ejObject to a dataframe with one row per record

#' @param x an ejObject
#' @param row.names NULL or a character vector giving the row names for the data frame. Missing values are not allowed.
#' @param optional not used.
#' @param ... other parameters passed to \code{\link{data.frame}}.
#'
#'
#' @return dataframe
#'
#' @method as.data.frame ejObject
#' @export
as.data.frame.ejObject <- function(x, row.names = NULL, optional = FALSE, ...){
	getAttList <- function(x){
		attList <- lapply(x, "[[", i="value")
		if(length(attList) > 0)
			names(attList) <- lapply(x, "[[", i="name")
		return(attList)
	}
	#metadata is added to the dataframe as part of the attrbutes
	metadata <- getAttList(x$metadata)
		
	rowList <- lapply(x$records, function(record){
				#grab the attributes as a row
				attList <- getAttList(record$attributes)
				
				#add the id
				recordElements <- c(list(id=record$id), attList)
				
				# now process the events
				#first we need unique names for each event
				uniqueEventNames <- make.names(sapply(record$events, "[[", i="name"), unique=TRUE)

				eventList <- lapply(1:length(record$events), function(i){
							event <- record$events[[i]]
							#grab the attributes as a row
							attList <- getAttList(event$attributes)
														
							#add the date and/or location
							placementList <- list()
							if(!is.null(event$date)){
								placementList[[paste(uniqueEventNames[[i]],"date",sep="_")]] <- event$date
							}
							if(class(event$location) == "SpatialPoints"){
								placementList[[paste(uniqueEventNames[[i]],"locationX",sep="_")]] <- event$location$x
								placementList[[paste(uniqueEventNames[[i]],"locationY",sep="_")]] <- event$location$y
								placementList[[paste(uniqueEventNames[[i]],"locationCRS",sep="_")]] <- proj4string(event$location)
							}
							#TODO: Should add polygons and lines here (even as GeoJSON txt column)
							if(length(attList)>0){
								return(data.frame(placementList, attList))
							} else {
								return(data.frame(placementList))
							}
						})	
			eventElements <- do.call(cbind,eventList)
			data.frame(recordElements,eventElements)	
			})
	#now bind up the results
	dF <- plyr::rbind.fill(rowList)
	
	#this is only included because the generic function has a row.names arg
	#and this is needed to pass check
	if(!missing(row.names)) {
		row.names(dF) <- row.names
	}
	
	#add in the metadata
	attributes(dF)$ejMetadata <- metadata
	
	#return the dataframe
	return(dF)
}

