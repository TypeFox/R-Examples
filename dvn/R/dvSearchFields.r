dvSearchFields <- function(dv = getOption('dvn'), browser=FALSE, ...){
	xml <- dvQuery(verb = "metadataSearchFields", query = NULL, dv = dv, browser=browser, ...)
	if(is.null(xml))
		invisible(NULL)
	if(browser==FALSE){
		searchterms <- xpathApply(xmlParse(xml),"//SearchableField")
		d <- data.frame(matrix(nrow=length(searchterms),ncol=2))
		names(d) <- c("fieldName","fieldDescription")
		for(i in 1:length(searchterms)){
			d$fieldName[i] <- xmlValue(xmlChildren(searchterms[[i]])$fieldName)
			d$fieldDescription[i] <- xmlValue(xmlChildren(searchterms[[i]])$fieldDescription)
		}
		return(d)
	}
}
