dvStudyStatement <-
function(   objectid, dv=getOption('dvn'), user=getOption('dvn.user'),
            pwd=getOption('dvn.pwd'), browser=FALSE, ...){
    if(is.null(user) | is.null(pwd))
        stop('Must specify username (`user`) and password (`pwd`)')
    xml <- dvDepositQuery(query=paste('statement/study/',objectid,sep=''), user=user, pwd=pwd, dv=dv, browser=browser, ...)
    if(is.null(xml))
		invisible(NULL)
	if(browser==FALSE){
        xmlout <- list()
        xml.list <- xmlToList(xml)
        xmlout$objectId <- objectid
        xmlout$id <- xml.list$id
        xmlout$title <- xml.list$title$text
        xmlout$author <- xml.list$author$name
        xmlout$updated <- xml.list$updated
        tmp <- sapply(xml.list[names(xml.list)=='category'], function(i) c(i[['.attrs']]['term'],i[['text']]))
        for(i in 1:ncol(tmp))
            xmlout[tmp[1,i]] <- tmp[2,i]
        rm(tmp)
        tmp <- xml.list[names(xml.list)=='entry']
        if(length(tmp)>1){
            resources <- t(sapply(tmp, function(i) c(i$content, title=i$title$text, summary=i$summary$text, updated=i$updated)))
            rownames(resources) <- seq(1,nrow(resources))
            rm(tmp)
            xmlout$files <- as.data.frame(resources, stringsAsFactors=FALSE)
            xmlout$files$fileId <- sapply(xmlout$files$src, function(i) strsplit(strsplit(i,'file/')[[1]][2],'/')[[1]][1])
        } else
            xmlout$files <- NULL
        xmlout$xml <- xml
        class(xmlout) <- c(class(xmlout),'dvStudyStatement')
        return(xmlout)
    }
}

print.dvStudyStatement <- function(x,...){
    cat('Study author: ',x$author,'\n')
    cat('Study title:  ',x$title,'\n')
    cat('ObjectId:     ',x$objectId,'\n')
    cat('Study URI:    ',x$id,'\n')
    cat('Last updated: ',x$updated,'\n')
    cat('Status:       ',x$latestVersionState,'\n')
    cat('Locked?       ',x$locked,'\n')
    if(!is.null(x$files)){
        cat('Files:\n')
        print(x$files[,c('src','type','updated','fileId')], right=FALSE)
    } else
        cat('Files:         None\n')
    invisible(x)
}
