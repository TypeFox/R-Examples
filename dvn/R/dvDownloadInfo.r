dvDownloadInfo <- function(fileid, dv = getOption('dvn'), browser=FALSE, ...){
	if(is.null(fileid))
		stop("Must specify 'fileid'")
	xml <- dvQuery(verb = "downloadInfo", query = fileid, dv = dv, browser=browser, ...)
	if(is.null(xml))
		invisible(NULL)
	else if(browser==FALSE){
		parsed <- xmlParse(xml)
        details <- list()
		services <- xpathApply(parsed,"//accessService")
		attrs <- xmlAttrs(xpathApply(parsed,"//studyFile")[[1]])
		details$fileId <- as.character(attrs[names(attrs)=="fileId"])
		details$fileName <- xmlValue(xpathApply(parsed,"//fileName")[[1]])
		details$fileMimeType <- xmlValue(xpathApply(parsed,"//fileMimeType")[[1]])
		details$fileSize <- xmlValue(xpathApply(parsed,"//fileSize")[[1]])
		details$authMethod <- xmlValue(xpathApply(parsed,"//authMethod")[[1]])
		attrs <- xmlAttrs(xpathApply(parsed,"//Authorization")[[1]])
		details$directAccess <- as.character(attrs[names(attrs)=="directAccess"])
		x <- xpathApply(parsed,"//accessPermissions")
		if(length(x)>0)
			details$accessPermissions <- xmlValue(x[[1]])
		x <- xpathApply(parsed,"//accessRestrictions")
		if(length(x)>0)
			details$accessRestrictions <- xmlValue(x[[1]])
		else
			details$accessRestrictions <- ""
		
		details$accessServicesSupported <- data.frame(matrix(nrow=length(services),ncol=4))
		names(details$accessServicesSupported) <- c("serviceName","serviceArgs","contentType","serviceDesc")
		for(i in 1:length(services)){
			tmp <- xmlChildren(services[[i]])
			details$accessServicesSupported$serviceName[i] <- xmlValue(tmp$serviceName)
			details$accessServicesSupported$serviceArgs[i] <- xmlValue(tmp$serviceArgs)
			details$accessServicesSupported$contentType[i] <- xmlValue(tmp$contentType)
			details$accessServicesSupported$serviceDesc[i] <- xmlValue(tmp$serviceDesc)
		}
        details$xml <- xml
        class(details) <- c(class(details),'dvDownloadInfo')
		return(details)
	}
}

print.dvDownloadInfo <- function(x,...){
    cat('File Name:      ',x$fileName,'\n')
    cat('File ID:        ',x$fileId,'\n')
    cat('File Type:      ',x$fileMimeType,'\n')
    cat('File Size:      ',x$fileSize,'\n')
    cat('Authentication: ',x$authMethod,'\n')
    cat('Direct Access?  ',x$directAccess,'\n',x$accessRestrictions,'\n')
    cat('Access Services:\n')
    print(x$accessServicesSupported)
    invisible(x)
}
