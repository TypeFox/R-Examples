NULL

#' @title Coerces an object to a \code{smet-class}  object
#' 
#' @description The method \code{as.smet} coerces an object or a charachter string to a SMET. If the object is 
#'
#' 
#' 
#' @param object the object to be coerced
#' @param mult,offset numeric vectors of unit multiplier and offset respectivaly
#' @param date.field   field name used for date and time. Default is \code{"timestamp"}, as used for \code{SMET} format. 
#' @param station.field field name used for station ID. Default is \code{"station_id"}, as used for \code{SMET} format. 
#' @param header.fields names used for the SMET header. Defaults are \code{c("longitude","latitude","station_id" ,"altitude","location")}
#' @param variables (optional) selection of variables hich can be exported to SMET formats. It is used only in case of two or more stations. 
#' @param metaparam metedata optional data frame containig meta info on variables. It can be entered as an attribute of \code{object}. See the structure of \code{metaparam} of \code{\link{meteofrance}}.It must contains \code{SMET_ID},\code{SMET_UNIT_MULTIPLIER},\code{SMET_UNIT_OFFSET} columns/fields.
#' @param force.multistation logical value. If it is \code{TRUE} the method is forced to return a list of SMET objects even in case of only one station 
#' @param file full filename of the reference SMET filename (not considered when \code{object} is \code{character} .
#' 
#' 
#' @param ... further arguments
#' 
#' @rdname as.smet
#' @export
#' 
#' @examples 
#' data(meteofrance)
#' 
#' 
#' ## Choose a particular station 
#' station_id <-  unique(meteofrance$station_id)[3]
#'
#' 
#' 
#'
# variables <-  c("latitude","longitude","station_id","altitude",
#     "location",
#' variables <- c("timestamp","DW","VW","TA","TD","RH","MFR_rr24",
#'       "MFR_tn12","MFR_tn24","MFR_tx12","MFR_tx24","HS","HS_fresh")
#' header <- c("longitude","latitude","station_id" ,"altitude","location")
#' names(header) <- header
#' 
############ attr(data,"header") <- lapply(X=header,FUN=function(x,data) {data[1,x]},data=data)
########## data <- data[,!(names(data) %in% header)]
##' 
#' 
#' data <- meteofrance[meteofrance$station_id==station_id,c(header,variables)]
#' metaparam <- attr(meteofrance,"metaparam")
#' metaparam <- metaparam[metaparam$SMET_ID %in% names(data),]

#' header <- lapply(X=header,FUN=function(x,data) {data[1,x]},data=data)
#' data <- data[,variables]
#' attr(data,"header") <- header
#' attr(data,"metaparam") <- metaparam
#' 
#' sm <- as.smet(data)
#' 
#' # In case of multiple station, it return a list of SMET-class objects: 
#' 
#' sm_multi <- as.smet(meteofrance,variables=variables)
#' 
#' 
#'   
#' 
#' 


as.smet <- function (object=NULL,...)  {
	
	
	return(standardGeneric("as.smet"))
	
}


NULL
#' 
#' @title as.smet
#' @description as.smet
#' @rdname as.smet
#' @method as.smet default
#' @aliases as.smet 
#' @export


setGeneric("as.smet",function (object,...)  {
	
	out <- NA
	warning("Object cannot be coerced as 'smet'") 
	return(out)
	
	
	
})


NULL
#'
#' @title as.smet
#' @description as.smet
#' @rdname as.smet
#' @method as.smet character
#' @aliases as.smet 
#' @export

setMethod("as.smet","character",function(object,...) {
	
	if (object %in% c("test","example")) {
		
		
		object <- system.file("examples/test.smet",package="RSMET")
	}		
			
			
	if (file.exists(object)==TRUE) {
		
	##	out <- smet(file=object,...)
		out <- RSMET::smet(file=object,...)
		
	}	else {
		
		value <- get(object)
		out <- RSMET::as.smet(value,...)
		
	}
			
	
	return(out)
	
	
})


NULL
#'
#' @title as.smet
#' @description as.smet
#' 
#' 
#' @rdname as.smet
#' @method as.smet data.frame
#' @aliases as.smet 
#' @export
#' 
#' 

setMethod("as.smet","data.frame",function(object,mult=NA,offset=NA,date.field="timestamp",station.field="station_id",header.fields=c("longitude","latitude","station_id" ,"altitude","location"),variables=NULL,force.multistation=FALSE,
				metaparam=attr(object,"metaparam"),file=NA,...) {
		
	### IN CASE VARIABLES AND HEADER FIELDS ARE SET AS A UNIQUE STRING!!		
	variables <- unlist(str_split(variables,","))
    header.fields <- unlist(str_split(header.fields,","))	
	############ MULTISTATION OPTION 
	multistation <- FALSE
	if (is.null(station.field)) station.field <- NA
	if (!is.na(station.field)) {
	
		if (station.field %in% names(object)) {
			station.names <- unique(object[,station.field])
			if (length(station.names)>1) { 
				multistation <- TRUE
			}
		}
	}	
	if (force.multistation==TRUE) multistation <- TRUE
	if (multistation==TRUE) {	
		
		
		header.fields <- union(header.fields,station.field)
		header.fields <- header.fields[header.fields %in% names(object)]
		names(header.fields) <- header.fields
		
		
		if (!is.null(variables)) {
			
			variables <- variables[variables %in% names(object)]
		##	print(c(header.fields,variables))
			object <- object[,c(header.fields,variables)]
			
		} else {
			
			
			variables <- names(object)[!(names(object) %in% header.fields)]
		}
		object <- split(object,object[,station.field])
		object <- base::lapply(X=object,FUN=function(x,header.fields,variables) {
					
				out <- x[,variables]	
				attr(out,"header") <- base::lapply(X=header.fields,FUN=function(i,x){x[1,i]},x=x)	
				
				
				return(out)
				
			},header.fields=header.fields,variables=variables)	
		
		out <- base::lapply(X=object,FUN=as.smet,mult=mult,offset=offset,date.field=date.field,station.field=station.field,header.fields=header.fields,metaparam=metaparam,file=NA,...)
		
		return(out)
		
	}
	
	
	############   END MULTISTATION OPTION  20151005 
	
	
	signature <- attr(object,"signature")
	header <- attr(object,"header")
   
	newheader <- list(...)$header
    newsignature <- list(...)$signature
	
	if (is.null(header)) {
		
		header <- newheader
		
	} else{
		
		###union_n <- union(names(header),names(newheader))
		
		header[names(newheader)] <- newheader
		
		
	}
	
	if (is.null(signature)) signature <- newsignature
	
	out <- as.smet("test")	
	
	
	out@header[names(header)] <- header 
	
	
	
	
   
	
   if (!is.null(signature)) out@signature <- signature
   
   object[is.na(object)] <- out@header$nodata
   
   
   
   
   out@data <- object
   
   out@header$fields <- names(out@data)
   
   out@header$units_offset <- array(0,ncol(object))
   out@header$units_multiplier <- array(1,ncol(object))
   
   
   names(out@header$fields) <- out@header$fields
   names(out@header$units_multiplier) <- out@header$fields
   names(out@header$units_offset) <- out@header$fields
   
   
   
  ##### NOT CORRECT!!! 
   
   if (!is.na(mult)) {
	   
	#   out@data[,names(mult)]  <- t(apply(X= out@data[,names(mult)],FUN=function(x,mult) {x/mult},mult=mult,MARGIN=1))
	   
	   out@header$units_multpier[names(mult)] <- mult
   }  
   
   if (!is.na(offset)) {
	   
	 #  out@data[,names(offset)]  <- t(apply(X= out@data[,names(offset)],FUN=function(x,mult) {x-offset},offset=offset,MARGIN=1))
	   
	   out@header$units_offset[names(offset)] <- offset
   } 
	
   #print("metaparam:")
   #str(metaparam)
   
   if (!is.null(metaparam)) {
	   
	   ids <- which(metaparam$SMET_ID %in% out@header$fields) 
	   out@header$units_multpier[metaparam$SMET_ID[ids]] <- metaparam$SMET_UNIT_MULTIPLIER[ids]
	   out@header$units_offset[metaparam$SMET_ID[ids]] <- metaparam$SMET_UNIT_OFFSET[ids]
	   
	   
	   
   }
   
   tt <- out@data[1,date.field]
   tzv <- round(as.numeric(as.POSIXlt(as.character(tt),tz="GMT")-tt,units="hours")/0.5)*0.5
   
   if (tzv>=0) {
	   
	   mntz <- round((tzv-trunc(tzv))*60)
	  	
	   if (mntz!=0) {
		   
		   tz <- sprintf("+%02d:%02d",trunc(tzv),mntz)
	   
	   } else {
		   
		   tz <- sprintf("+%02d",trunc(tzv))
		   
	   }
	   
	   
	   
	   
   } else {
	   
	   mntz <- abs(round((tzv-trunc(tzv))*60))
	   
	   if (mntz!=0) {
		   
		   tz <- sprintf("-%02d:%02d",-trunc(tzv),mntz)
		   
	   } else {
		   
		   tz <- sprintf("-%02d",-trunc(tzv))
		   
	   }
	   
   }
  
   out@header$tz <- tz
   
   out@file <- as.character(file)
	
	#### SORTING DATA  ### EC 20151104
	
	if (date.field %in% names(out@data)) {
		
		
		row.names(out@data) <- sprintf("%09d",1:nrow(out@data))
		times <- out@data[,date.field]
		names(times) <- row.names(out@data)
		
		out@data <- out@data[names(sort(times)),]
		
		
		
		
		
		
	}
	
	
	#####
	
	
	return(out)
	
	
})


NULL
#'
#' @title as.smet
#' @description as.smet
#' 
#' 
#' @rdname as.smet
#' @method as.smet list
#' @aliases as.smet 
#' @export
#' 
#' 


setMethod("as.smet","list",function(object,...) {base::lapply(X=object,FUN=RSMET::as.smet,...)})
					   


	
	
NULL
#'
#' @title as.smet
#' @description as.smet
#' 
#' 
#' @rdname as.smet
#' @method as.smet smet
#' @aliases as.smet 
#' @export
#' 
#' 
	
	
	setMethod("as.smet","smet",function(object,...) { 
				
				
				args <- base::list(...)
				slotnames <- names(getSlots("smet"))
				
				slotnames <- slotnames[slotnames %in% names(args)]
				
				for (it in slotnames) {
					
					slot(object,it) <- args[[it]]
				}
				
				
				return(object)})
	
	




































