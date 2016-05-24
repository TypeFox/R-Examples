dvTermsOfUse <- function(xml){
    if(inherits(xml,'dvServiceDoc')){
        out <- lapply(xml$dataverses$collectionPolicy, function(x){
            writeLines(x, tmp <- tempfile(fileext='.html'))
            browseURL(tmp)
            return(tmp)
        })
        Sys.sleep(1)
        unlink(out)
        return(invisible(xml))
    }
    else if(inherits(xml,'dvMetadata')){
        f <- attributes(xml)$formatName
        if(is.null(f) || !f=='ddi')
            stop("Object of class 'dvMetadata' must have format 'ddi'.")
        tmp <- xmlToList(xml)$stdyDscr$dataAccs$notes$text
        writeLines(tmp, out <- tempfile(fileext='.html'))
        browseURL(out)
        Sys.sleep(1)
        unlink(out)
        return(invisible(xml))
    }
    else
        stop('Currently only objects of class `dvServiceDoc` are supported.')
}
