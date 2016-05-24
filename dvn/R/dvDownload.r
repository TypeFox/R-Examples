dvDownload <- function(fileid, query=NULL, dv = getOption('dvn'), browser=FALSE, ...){
	if(is.null(fileid))
		stop("Must specify 'fileid'")
	direct <- dvDownloadInfo(fileid, dv = dv, ...)
	if(is.null(direct))
		stop("downloadInfo unavailable")
	if(direct$directAccess=="false")
		stop(direct$accessRestrictions,"\nData cannot be accessed directly...try browsing URI from dvExtractFileIds(dvMetadata())")
	if(is.null(query)){
		xml <- dvQuery(verb = "download", query = fileid, dv = dv, browser=browser)
		return(xml)
	}
	else{
		if(is.list(query)){
			qtmp <- ""
			for(i in 1:length(query)){
				qtmp <- paste(qtmp,names(query)[i],"=",query[[i]],sep="")
				if(i<length(query))
					qtmp <- paste(qtmp,"&",sep="")
			}
			query <- qtmp
		}
		else if(is.character(query))
			query <- query
		else
			stop("Must specify query as named list or character string")
		xml <- dvQuery(verb = "download", query = paste(fileid,query,sep=""), dv = dv, browser=browser, ...)
        if(is.null(xml))
            invisible(NULL)
        else if(browser==FALSE)
            return(xml)
    }
}
