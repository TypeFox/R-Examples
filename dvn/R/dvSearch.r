dvSearch <- function(query, boolean='AND', dv = getOption('dvn'), browser=FALSE, ...){
	if(is.null(query))
		stop("Must specify query as named list or character string")
	else{
		if(is.list(query)){
            if(length(query)==1 && is.null(names(query)))
			    query <- query
            else{
                qtmp <- ""
                for(i in 1:length(query)){
                    qtmp <- paste(qtmp,names(query)[i],":",curlEscape(query[[i]]),sep="")
                    if(i<length(query))
                        qtmp <- paste(qtmp,"%20",boolean,"%20",sep="")
                }
                query <- qtmp
            }
		}
		else if(is.character(query))
			query <- paste(dvSearchFields()$fieldName,curlEscape(query[1]),sep=':', collapse='%20OR%20')
		else
			stop("Must specify query as named list or character string")
	}
	xml <- dvQuery(verb = "metadataSearch", query = query, dv = dv, browser=browser, ...)
	if(is.null(xml))
		return(NULL)
	else if(browser==FALSE){
		results <- unlist(xpathApply(xmlParse(xml),"//study", fun=xmlAttrs))
		if(length(results))
			d <- data.frame(objectId=results)
		else
			d <- NULL
		message(if(is.null(d)) '0' else nrow(d), ' search results returned\n')
		return(d)
	}
}
