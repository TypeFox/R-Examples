dvMetadata <- function(objectid, format.type=NULL, include=NULL, exclude=NULL,
						dv = getOption('dvn'), browser=FALSE, ...){
	if(is.null(format.type)){
		query <- objectid
		if(!is.null(include))
			query <- paste(query,"/?partialInclude=",include,sep="")
		if(!is.null(exclude))
			query <- paste(query,"/?partialExclude=",exclude,sep="")
	}
	else{
		query <- paste(objectid,"/?formatType=",format.type,sep="")
		if(!is.null(include))
			query <- paste(query,"&partialInclude=",include,sep="")
		if(!is.null(exclude))
			query <- paste(query,"&partialExclude=",exclude,sep="")
	}
	xml <- dvQuery(verb = "metadata", query = query, dv = dv, browser=browser, ...)
	if(is.null(xml))
		invisible(NULL)
	else if(browser==FALSE){
        class(xml) <- c(class(xml),'dvMetadata')
        attr(xml,'formatName') <- ifelse(!is.null(format.type),format.type,'ddi')
		return(xml)
    }
}
