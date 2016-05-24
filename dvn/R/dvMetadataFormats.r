dvMetadataFormats <- function(objectid, dv = getOption('dvn'), browser=FALSE, ...){
	xml <- dvQuery(verb = "metadataFormatsAvailable", query = objectid, dv = dv, browser=browser, ...)
	if(is.null(xml))
		invisible(NULL)
	else if(browser==FALSE){
		searchterms <- xpathApply(xmlParse(xml),"//formatAvailable")
		if(length(searchterms)>0){
			d <- data.frame(matrix(nrow=length(searchterms),ncol=3))
			names(d) <- c("formatName","formatSchema","formatMime")
			for(i in 1:length(searchterms)){
				d$formatName[i] <- xmlValue(xmlChildren(searchterms[[i]])$formatName)
				d$formatSchema[i] <- xmlValue(xmlChildren(searchterms[[i]])$formatSchema)
				d$formatMime[i] <- xmlValue(xmlChildren(searchterms[[i]])$formatMime)
			}
			return(d)
		}
		else
			return(NULL)
	}
}
