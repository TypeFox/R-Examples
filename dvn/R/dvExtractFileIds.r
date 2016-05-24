dvExtractFileIds <- function(xml){
	if(is.null(xmlChildren(xmlParse(xml))$codeBook))
        stop("Metadata format not currently supported. Must be 'ddi'.")
    nodes <- xmlChildren(xmlChildren(xmlParse(xml))$codeBook)
	
    # data files
    dscrs <- nodes[names(nodes)=="fileDscr"]
	d <- 
    lapply(dscrs, function(i){
        attrs <- xmlAttrs(i)
        fileName <- xmlValue(xmlChildren(xmlChildren(i)$fileTxt)$fileName)
        fileId <- attrs[names(attrs)=="ID"]
		fileId <- substring(fileId,2,nchar(fileId)) # remove leading 'f' from fileId
		URI <- attrs[names(attrs)=="URI"]
		dims <- xmlChildren(xmlChildren(xmlChildren(i)$fileTxt)$dimensns)
		caseQnty <- xmlValue(dims$caseQnty)
		varQnty <- xmlValue(dims$varQnty)
        notes <- xmlChildren(i)[names(xmlChildren(i))=='notes']
        fingerprint <- 
            xmlValue(notes[sapply(notes, function(i)
                'Universal Numeric Fingerprint' %in% xmlAttrs(i))]$notes)
        return(list(fileName = fileName,
                    fileId = fileId,
                    UNF = fingerprint,
                    URI = URI,
                    caseQnty = caseQnty,
                    varQnty = varQnty))
    })
    d <- do.call(rbind.data.frame, d)
    rownames(d) <- NULL
    
    # other files
    other <- nodes[names(nodes)=="otherMat"]
    e <- 
    lapply(other, function(i){
        URI <- xmlAttrs(i)['URI']
        fileName <- xmlValue(xmlChildren(i)$labl)
        fileId <- strsplit(URI, 'fileId=')[[1]][2]
        return(list(fileName = fileName,
                    fileId = fileId,
                    URI = URI))
    })
    e <- do.call(rbind.data.frame, e)
    rownames(e) <- NULL
    
    if(length(d)>0 && length(e)>0)
        return(merge(d, e, all=TRUE)[,c('fileName','fileId','UNF','caseQnty','varQnty','URI')])
    else if(length(d)>0)
        return(d)
    else
        return(e)
}
